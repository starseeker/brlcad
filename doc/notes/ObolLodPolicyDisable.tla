----------------------- MODULE ObolLodPolicyDisable -----------------------
\* Focused ownership contract for toggling automatic mesh LoD.
\*
\* A policy-off transition retires every automatic producer, cursor, and
\* capacity-labelled frame atomically.  It does not delete the immutable
\* presentation already on screen, cancel an independent database provider,
\* or suppress the one presentation-only repaint needed to show the policy
\* change.  The model deliberately abstracts geometry and numeric quality.

EXTENDS Naturals, TLC

CONSTANT MaxPolicyEpoch

WorkerStates == {"idle", "queued", "inflight", "result"}

VARIABLES automatic,
          policyEpoch,
          retainedPresentation,
          workerState,
          submissionPending,
          allocationPending,
          presentationPending,
          interactionActive,
          compactionPending,
          capacityRenderPending,
          repaintPending,
          providerPending

vars == <<automatic, policyEpoch, retainedPresentation, workerState,
          submissionPending, allocationPending, presentationPending,
          interactionActive, compactionPending, capacityRenderPending,
          repaintPending, providerPending>>

TypeOK ==
    /\ automatic \in BOOLEAN
    /\ policyEpoch \in 0..MaxPolicyEpoch
    /\ retainedPresentation \in BOOLEAN
    /\ workerState \in WorkerStates
    /\ submissionPending \in BOOLEAN
    /\ allocationPending \in BOOLEAN
    /\ presentationPending \in BOOLEAN
    /\ interactionActive \in BOOLEAN
    /\ compactionPending \in BOOLEAN
    /\ capacityRenderPending \in BOOLEAN
    /\ repaintPending \in BOOLEAN
    /\ providerPending \in BOOLEAN

Init ==
    /\ automatic = TRUE
    /\ policyEpoch = 0
    /\ retainedPresentation = TRUE
    /\ workerState \in (WorkerStates \ {"idle"})
    /\ submissionPending \in BOOLEAN
    /\ allocationPending \in BOOLEAN
    /\ presentationPending \in BOOLEAN
    /\ interactionActive \in BOOLEAN
    /\ compactionPending \in BOOLEAN
    /\ capacityRenderPending \in BOOLEAN
    /\ repaintPending \in BOOLEAN
    /\ providerPending \in BOOLEAN

\* Odd transitions disable; even transitions re-enable.  Re-enable starts a
\* fresh worker/submission epoch but continues to display the retained image.
TogglePolicy ==
    /\ policyEpoch < MaxPolicyEpoch
    /\ policyEpoch' = policyEpoch + 1
    /\ automatic' = ~automatic
    /\ retainedPresentation' = retainedPresentation
    /\ repaintPending' = TRUE
    /\ providerPending' = providerPending
    /\ IF automatic
          THEN /\ workerState' = "idle"
               /\ submissionPending' = FALSE
               /\ allocationPending' = FALSE
               /\ presentationPending' = FALSE
               /\ interactionActive' = FALSE
               /\ compactionPending' = FALSE
               /\ capacityRenderPending' = FALSE
          ELSE /\ workerState' = "queued"
               /\ submissionPending' = TRUE
               /\ allocationPending' = FALSE
               /\ presentationPending' = FALSE
               /\ interactionActive' = FALSE
               /\ compactionPending' = FALSE
               /\ capacityRenderPending' = TRUE

\* Automatic policy may accumulate any combination of bounded work before a
\* later toggle.  The disable transition must retire all combinations, not
\* only the state produced by one expected ordering.
ArmAutomaticWork ==
    /\ automatic
    /\ \/ /\ workerState = "idle"
           /\ workerState' \in (WorkerStates \ {"idle"})
           /\ UNCHANGED <<submissionPending, allocationPending,
                           presentationPending, interactionActive,
                           compactionPending, capacityRenderPending>>
       \/ /\ submissionPending' = TRUE
           /\ UNCHANGED <<workerState, allocationPending,
                           presentationPending, interactionActive,
                           compactionPending, capacityRenderPending>>
       \/ /\ allocationPending' = TRUE
           /\ UNCHANGED <<workerState, submissionPending,
                           presentationPending, interactionActive,
                           compactionPending, capacityRenderPending>>
       \/ /\ presentationPending' = TRUE
           /\ UNCHANGED <<workerState, submissionPending,
                           allocationPending, interactionActive,
                           compactionPending, capacityRenderPending>>
       \/ /\ interactionActive' = TRUE
           /\ UNCHANGED <<workerState, submissionPending,
                           allocationPending, presentationPending,
                           compactionPending, capacityRenderPending>>
       \/ /\ compactionPending' = TRUE
           /\ UNCHANGED <<workerState, submissionPending,
                           allocationPending, presentationPending,
                           interactionActive, capacityRenderPending>>
       \/ /\ capacityRenderPending' = TRUE
           /\ UNCHANGED <<workerState, submissionPending,
                           allocationPending, presentationPending,
                           interactionActive, compactionPending>>
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    repaintPending, providerPending>>

CompleteProvider ==
    /\ providerPending
    /\ providerPending' = FALSE
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    workerState, submissionPending, allocationPending,
                    presentationPending, interactionActive,
                    compactionPending, capacityRenderPending, repaintPending>>

PresentRepaint ==
    /\ repaintPending
    /\ repaintPending' = FALSE
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    workerState, submissionPending, allocationPending,
                    presentationPending, interactionActive,
                    compactionPending, capacityRenderPending, providerPending>>

Next == TogglePolicy \/ ArmAutomaticWork \/ CompleteProvider \/ PresentRepaint

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(TogglePolicy)
    /\ WF_vars(CompleteProvider)
    /\ WF_vars(PresentRepaint)

NoAutomaticOwner ==
    /\ workerState = "idle"
    /\ ~submissionPending
    /\ ~allocationPending
    /\ ~presentationPending
    /\ ~interactionActive
    /\ ~compactionPending
    /\ ~capacityRenderPending

RetainedPresentationSurvives == retainedPresentation

DisabledOwnsNoAutomaticWork == ~automatic => NoAutomaticOwner

DisabledRenderIsPresentationOnly == ~automatic => ~capacityRenderPending

EventuallyQuiescentAfterFinalDisable ==
    <> (policyEpoch = MaxPolicyEpoch /\ ~automatic /\ NoAutomaticOwner /\
        ~providerPending /\ ~repaintPending)

=============================================================================
