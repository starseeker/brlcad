-------------------- MODULE ObolActiveProducerDemand --------------------
\* One immutable asset producer may outlive several view/policy demands.
\* This model checks only that ownership seam: current-demand results may be
\* applied, results overtaken by a demand or compact-population generation are
\* superseded rather than provider failures, and an exact-demand failure cannot
\* survive a demand epoch change.  Consuming a superseded result must leave
\* NeedCurrentDemand true so fair submission installs its successor witness.

EXTENDS Naturals, TLC

CONSTANT MaxDemand

Demands == 0..MaxDemand
NoneDemand == MaxDemand + 1
OptionalDemand == 0..NoneDemand
ProducerStates == {"idle", "active"}
ResultStates == {"none", "ready", "superseded", "failure"}

VARIABLES inputOpen,
          demand,
          producer,
          providerDemand,
          resultState,
          resultDemand,
          previewDemand,
          presentedDemand,
          failureDemand

vars == <<inputOpen, demand, producer, providerDemand, resultState,
          resultDemand, previewDemand, presentedDemand, failureDemand>>

TypeOK ==
    /\ inputOpen \in BOOLEAN
    /\ demand \in Demands
    /\ producer \in ProducerStates
    /\ providerDemand \in Demands
    /\ resultState \in ResultStates
    /\ resultDemand \in Demands
    /\ previewDemand \in OptionalDemand
    /\ presentedDemand \in OptionalDemand
    /\ failureDemand \in OptionalDemand

Init ==
    /\ inputOpen = TRUE
    /\ demand = 0
    /\ producer = "idle"
    /\ providerDemand = 0
    /\ resultState = "none"
    /\ resultDemand = 0
    /\ previewDemand = NoneDemand
    /\ presentedDemand = NoneDemand
    /\ failureDemand = NoneDemand

MoveDemand ==
    /\ inputOpen
    /\ demand < MaxDemand
    /\ demand' = demand + 1
    /\ failureDemand' = NoneDemand
    /\ UNCHANGED <<inputOpen, producer, providerDemand, resultState,
                    resultDemand, previewDemand, presentedDemand>>

CloseInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ UNCHANGED <<demand, producer, providerDemand, resultState,
                    resultDemand, previewDemand, presentedDemand,
                    failureDemand>>

NeedCurrentDemand ==
    /\ presentedDemand # demand
    /\ failureDemand # demand

Submit ==
    /\ producer = "idle"
    /\ resultState = "none"
    /\ NeedCurrentDemand
    /\ producer' = "active"
    /\ providerDemand' = demand
    /\ UNCHANGED <<inputOpen, demand, resultState, resultDemand,
                    previewDemand, presentedDemand, failureDemand>>

\* Intermediate immutable geometry is rebound immediately before publication.
\* A later demand may overtake the queued preview; the owner then rejects it.
PublishPreview ==
    /\ producer = "active"
    /\ previewDemand' = demand
    /\ UNCHANGED <<inputOpen, demand, producer, providerDemand,
                    resultState, resultDemand, presentedDemand,
                    failureDemand>>

\* The provider refreshes view-dependent page selection after stable asset
\* construction.  This does not change asset identity or restart the worker.
RefreshDemand ==
    /\ producer = "active"
    /\ providerDemand # demand
    /\ providerDemand' = demand
    /\ UNCHANGED <<inputOpen, demand, producer, resultState, resultDemand,
                    previewDemand, presentedDemand, failureDemand>>

Complete ==
    /\ producer = "active"
    /\ producer' = "idle"
    /\ resultDemand' = providerDemand
    /\ resultState' = IF providerDemand = demand THEN "ready" ELSE "superseded"
    /\ UNCHANGED <<inputOpen, demand, providerDemand, previewDemand,
                    presentedDemand, failureDemand>>

ProviderFailure ==
    /\ producer = "active"
    /\ providerDemand = demand
    /\ producer' = "idle"
    /\ resultDemand' = demand
    /\ resultState' = "failure"
    /\ UNCHANGED <<inputOpen, demand, providerDemand, previewDemand,
                    presentedDemand, failureDemand>>

Consume ==
    /\ resultState # "none"
    /\ presentedDemand' =
        IF resultState = "ready" /\ resultDemand = demand
        THEN demand ELSE presentedDemand
    /\ failureDemand' =
        IF resultState = "failure" /\ resultDemand = demand
        THEN demand ELSE failureDemand
    /\ resultState' = "none"
    /\ UNCHANGED <<inputOpen, demand, producer, providerDemand,
                    resultDemand, previewDemand>>

Next ==
    \/ MoveDemand
    \/ CloseInput
    \/ Submit
    \/ PublishPreview
    \/ RefreshDemand
    \/ Complete
    \/ ProviderFailure
    \/ Consume

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(Submit)
    /\ WF_vars(Complete)
    /\ WF_vars(ProviderFailure)
    /\ WF_vars(Consume)

MonotonePublishedDemand ==
    /\ (previewDemand # NoneDemand => previewDemand <= demand)
    /\ (presentedDemand # NoneDemand => presentedDemand <= demand)

FailureIsCurrentDemand ==
    failureDemand = NoneDemand \/ failureDemand = demand

\* Service normalization and source result routing classify an overtaken demand
\* or retired compact population as superseded; only an actual provider failure
\* may create exact-demand failure state.
SupersededResultIsNotFailure ==
    resultState = "superseded" => failureDemand # resultDemand

Terminal ==
    /\ ~inputOpen
    /\ producer = "idle"
    /\ resultState = "none"
    /\ (presentedDemand = demand \/ failureDemand = demand)

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyTerminalAfterInputCloses ==
    [](~inputOpen => <>Terminal)

=============================================================================
