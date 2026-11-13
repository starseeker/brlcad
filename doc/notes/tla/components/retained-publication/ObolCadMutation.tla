------------------------- MODULE ObolCadMutation -------------------------
\* Atomic retained-scene mutation after side-effect-free validation.
\*
\* Validation and staging may happen outside the presentation node.  Every
\* nonempty typed mutation journal -- including a sparse cut-only journal --
\* enters one RAII update scope; the renderer-visible scene changes only when
\* that scope closes.  Empty, invalid, or resource-denied candidates cannot
\* open a scope or change the preceding scene.  Resource denial represents
\* failure while allocating a complete replacement candidate before the
\* no-throw publication swap.  notifiedScene is the scene last exposed by
\* the node notification/touch edge; it must advance in the same CloseCommit
\* as liveScene.  Transaction framing therefore cannot be part of swappable
\* retained payload.

EXTENDS Naturals, TLC

CONSTANT Scenes, MaxBatches

CandidateStates == {"none", "valid", "invalid", "unavailable"}
MutationKinds == {"replacement", "part", "instance", "style", "cut",
                  "part-removal", "instance-removal", "instance-set",
                  "draw-mode", "semantic"}

VARIABLES inputClosed,
          liveScene,
          notifiedScene,
          candidateScene,
          candidateState,
          candidateKinds,
          updateOpen,
          baselineScene,
          workingScene,
          remainingBatches

vars == <<inputClosed, liveScene, notifiedScene,
          candidateScene, candidateState,
          candidateKinds,
          updateOpen, baselineScene, workingScene, remainingBatches>>

TypeOK ==
    /\ inputClosed \in BOOLEAN
    /\ liveScene \in Scenes
    /\ notifiedScene \in Scenes
    /\ candidateScene \in Scenes \union {"none"}
    /\ candidateState \in CandidateStates
    /\ candidateKinds \subseteq MutationKinds
    /\ updateOpen \in BOOLEAN
    /\ baselineScene \in Scenes \union {"none"}
    /\ workingScene \in Scenes \union {"none"}
    /\ remainingBatches \in 0..MaxBatches

Init ==
    /\ inputClosed = FALSE
    /\ liveScene \in Scenes
    /\ notifiedScene = liveScene
    /\ candidateScene = "none"
    /\ candidateState = "none"
    /\ candidateKinds = {}
    /\ updateOpen = FALSE
    /\ baselineScene = "none"
    /\ workingScene = "none"
    /\ remainingBatches = 0

StageValid ==
    /\ ~inputClosed
    /\ candidateState = "none"
    /\ ~updateOpen
    /\ \E scene \in Scenes, kinds \in SUBSET MutationKinds:
        /\ kinds # {}
        /\ candidateScene' = scene
        /\ candidateState' = "valid"
        /\ candidateKinds' = kinds
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

StageInvalid ==
    /\ ~inputClosed
    /\ candidateState = "none"
    /\ ~updateOpen
    /\ \E scene \in Scenes, kinds \in SUBSET MutationKinds:
        /\ kinds # {}
        /\ candidateScene' = scene
        /\ candidateState' = "invalid"
        /\ candidateKinds' = kinds
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

StageUnavailable ==
    /\ ~inputClosed
    /\ candidateState = "none"
    /\ ~updateOpen
    /\ \E scene \in Scenes, kinds \in SUBSET MutationKinds:
        /\ kinds # {}
        /\ candidateScene' = scene
        /\ candidateState' = "unavailable"
        /\ candidateKinds' = kinds
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<liveScene, notifiedScene,
                    candidateScene, candidateState,
                    candidateKinds, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

RejectInvalid ==
    /\ candidateState = "invalid"
    /\ ~updateOpen
    /\ candidateScene' = "none"
    /\ candidateState' = "none"
    /\ candidateKinds' = {}
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

RejectUnavailable ==
    /\ candidateState = "unavailable"
    /\ ~updateOpen
    /\ candidateScene' = "none"
    /\ candidateState' = "none"
    /\ candidateKinds' = {}
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, updateOpen,
                    baselineScene, workingScene, remainingBatches>>

BeginCommit ==
    /\ candidateState = "valid"
    /\ candidateKinds # {}
    /\ ~updateOpen
    /\ updateOpen' = TRUE
    /\ baselineScene' = liveScene
    /\ workingScene' = candidateScene
    /\ remainingBatches' \in 1..MaxBatches
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, candidateScene,
                    candidateState, candidateKinds>>

ApplyBatch ==
    /\ updateOpen
    /\ remainingBatches > 0
    /\ remainingBatches' = remainingBatches - 1
    /\ UNCHANGED <<inputClosed, liveScene, notifiedScene, candidateScene,
                    candidateState, candidateKinds, updateOpen,
                    baselineScene, workingScene>>

CloseCommit ==
    /\ updateOpen
    /\ remainingBatches = 0
    /\ liveScene' = workingScene
    /\ notifiedScene' = workingScene
    /\ candidateScene' = "none"
    /\ candidateState' = "none"
    /\ candidateKinds' = {}
    /\ updateOpen' = FALSE
    /\ baselineScene' = "none"
    /\ workingScene' = "none"
    /\ remainingBatches' = 0
    /\ UNCHANGED inputClosed

Quiescent ==
    /\ inputClosed
    /\ candidateState = "none"
    /\ ~updateOpen
    /\ UNCHANGED vars

Next ==
    \/ StageValid
    \/ StageInvalid
    \/ StageUnavailable
    \/ CloseInput
    \/ RejectInvalid
    \/ RejectUnavailable
    \/ BeginCommit
    \/ ApplyBatch
    \/ CloseCommit
    \/ Quiescent

Spec == Init /\ [][Next]_vars
        /\ WF_vars(RejectInvalid)
        /\ WF_vars(RejectUnavailable)
        /\ WF_vars(BeginCommit)
        /\ WF_vars(ApplyBatch)
        /\ WF_vars(CloseCommit)

CandidateShape ==
    (candidateState = "none") =
        (candidateScene = "none" /\ candidateKinds = {})

OpenScopeHasValidCandidate ==
    updateOpen =>
        /\ candidateState = "valid"
        /\ candidateKinds # {}
        /\ candidateScene \in Scenes
        /\ baselineScene \in Scenes
        /\ workingScene = candidateScene

OpenScopeOwnsFiniteMechanicalWork ==
    updateOpen => remainingBatches \in 0..MaxBatches

LiveSceneIsStableInsideScope ==
    updateOpen => liveScene = baselineScene

PublishedSceneWasNotified ==
    liveScene = notifiedScene

UncommittableCandidateCannotMutate ==
    candidateState \in {"invalid", "unavailable"} =>
        /\ ~updateOpen
        /\ baselineScene = "none"
        /\ workingScene = "none"
        /\ remainingBatches = 0

NonemptyValidJournalHasCommitOwner ==
    candidateState = "valid" /\ ~updateOpen => ENABLED BeginCommit

EventuallyDrainedAfterStableInput ==
    [](inputClosed => <>(candidateState = "none" /\ ~updateOpen))

=============================================================================
