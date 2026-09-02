---------------- MODULE ObolControlLifecycleComposition ----------------
\* Composition contract for the controller lifecycle seams which are only
\* assumptions in the data-plane models: policy enablement, provider
\* registration versus provider work, view-demand supersession, host pumping,
\* exact-frame acceptance, and terminal publication.
\* A hard deadline miss is also a first-class lifecycle edge: the first miss
\* of a frozen control problem must create one quality-recovery owner and a
\* fresh exact-frame obligation.  It cannot merely discard the interrupted
\* frame and leave a nonterminal view with no producer.
\*
\* This is intentionally smaller than ObolProgressivePipeline.  Quality is one
\* finite owner and geometry is absent.  Its purpose is to reject a lifecycle
\* implementation in which a standing fact has no host witness, policy-off
\* camera bookkeeping creates automatic work, an idle registered provider is
\* treated as discovery, or a stale frame settles the current view.

EXTENDS Naturals, TLC

CONSTANTS MaxViewRevision, MaxPolicyChanges, MaxImportanceRequests,
          MaxSemanticChanges, MaxDemandInterruptions

ASSUME MaxViewRevision > 0
ASSUME MaxPolicyChanges > 0
ASSUME MaxImportanceRequests > 0
ASSUME MaxSemanticChanges > 0
ASSUME MaxDemandInterruptions > 0

RenderStates == {"idle", "pending", "inflight"}
Outcomes == {"active", "ready", "constrained"}

VARIABLES automatic,
          policyChanges,
          importanceRequests,
          semanticChanges,
          viewRevision,
          committedRevision,
          semanticRevision,
          committedSemanticRevision,
          providerRegistered,
          providersClosed,
          providerPending,
          demandPending,
          demandCursorActive,
          demandInterruptions,
          importanceCensusPending,
          qualityPending,
          exactFrameRequired,
          capacityFrameRequired,
          exactFrameQueued,
          pumpPending,
          renderState,
          renderTargetRevision,
          renderTargetSemanticRevision,
          frameSuperseded,
          constraintWitness,
          deadlineRecoveryPending,
          deadlineRetryAvailable,
          retainedPresentation

vars == <<automatic, policyChanges, importanceRequests, semanticChanges,
          viewRevision, committedRevision,
          semanticRevision, committedSemanticRevision,
          providerRegistered, providersClosed, providerPending,
          demandPending, demandCursorActive, demandInterruptions,
          importanceCensusPending, qualityPending,
          exactFrameRequired, capacityFrameRequired,
          exactFrameQueued, pumpPending, renderState,
          renderTargetRevision, renderTargetSemanticRevision,
          frameSuperseded, constraintWitness,
          deadlineRecoveryPending, deadlineRetryAvailable,
          retainedPresentation>>

PumpReasons ==
    providerPending \/ demandPending \/ importanceCensusPending \/
    qualityPending \/
    deadlineRecoveryPending \/
    (exactFrameRequired /\ ~exactFrameQueued)

ForegroundWork ==
    PumpReasons \/ renderState # "idle"

Outcome ==
    IF ForegroundWork \/ exactFrameRequired \/
       committedRevision # viewRevision \/
       committedSemanticRevision # semanticRevision
    THEN "active"
    ELSE IF constraintWitness THEN "constrained" ELSE "ready"

TypeOK ==
    /\ automatic \in BOOLEAN
    /\ policyChanges \in 0..MaxPolicyChanges
    /\ importanceRequests \in 0..MaxImportanceRequests
    /\ semanticChanges \in 0..MaxSemanticChanges
    /\ viewRevision \in 0..MaxViewRevision
    /\ committedRevision \in 0..MaxViewRevision
    /\ committedRevision <= viewRevision
    /\ semanticRevision \in 0..MaxSemanticChanges
    /\ committedSemanticRevision \in 0..MaxSemanticChanges
    /\ committedSemanticRevision <= semanticRevision
    /\ providerRegistered \in BOOLEAN
    /\ providersClosed \in BOOLEAN
    /\ providerPending \in BOOLEAN
    /\ providerPending => providerRegistered
    /\ providersClosed => ~providerRegistered
    /\ demandPending \in BOOLEAN
    /\ demandCursorActive \in BOOLEAN
    /\ demandInterruptions \in 0..MaxDemandInterruptions
    /\ demandCursorActive => automatic /\ demandPending
    /\ importanceCensusPending \in BOOLEAN
    /\ qualityPending \in BOOLEAN
    /\ ~(demandPending /\ qualityPending)
    /\ exactFrameRequired \in BOOLEAN
    /\ capacityFrameRequired \in BOOLEAN
    /\ capacityFrameRequired => exactFrameRequired
    /\ committedSemanticRevision # semanticRevision => exactFrameRequired
    /\ exactFrameQueued \in BOOLEAN
    /\ exactFrameQueued => exactFrameRequired
    /\ pumpPending \in BOOLEAN
    /\ renderState \in RenderStates
    /\ renderTargetRevision \in 0..MaxViewRevision
    /\ renderTargetSemanticRevision \in 0..MaxSemanticChanges
    /\ frameSuperseded \in BOOLEAN
    /\ constraintWitness \in BOOLEAN
    /\ deadlineRecoveryPending \in BOOLEAN
    /\ deadlineRetryAvailable \in BOOLEAN
    /\ deadlineRecoveryPending => qualityPending
    /\ retainedPresentation \in BOOLEAN
    /\ Outcome \in Outcomes

Init ==
    /\ automatic = TRUE
    /\ policyChanges = 0
    /\ importanceRequests = 0
    /\ semanticChanges = 0
    /\ viewRevision = 0
    /\ committedRevision = 0
    /\ semanticRevision = 0
    /\ committedSemanticRevision = 0
    /\ providerRegistered \in BOOLEAN
    /\ providersClosed = FALSE
    /\ providerPending = providerRegistered
    /\ demandPending = ~providerRegistered
    /\ demandCursorActive = FALSE
    /\ demandInterruptions = 0
    /\ importanceCensusPending = FALSE
    /\ qualityPending = FALSE
    /\ exactFrameRequired = TRUE
    /\ capacityFrameRequired = TRUE
    /\ exactFrameQueued = FALSE
    /\ pumpPending = TRUE
    /\ renderState = "idle"
    /\ renderTargetRevision = 0
    /\ renderTargetSemanticRevision = 0
    /\ frameSuperseded = FALSE
    /\ constraintWitness = FALSE
    /\ deadlineRecoveryPending = FALSE
    /\ deadlineRetryAvailable = TRUE
    /\ retainedPresentation = TRUE

RegisterProvider ==
    /\ ~providersClosed
    /\ ~providerRegistered
    /\ providerRegistered' = TRUE
    /\ providerPending' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providersClosed, demandPending,
                    demandCursorActive, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

PublishProviderWork ==
    /\ providerRegistered
    /\ ~providerPending
    /\ providerPending' = TRUE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    demandPending, demandCursorActive, demandInterruptions,
                    importanceCensusPending,
                    qualityPending, exactFrameRequired,
                    capacityFrameRequired, exactFrameQueued,
                    renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

CloseProviders ==
    /\ ~providersClosed
    /\ providersClosed' = TRUE
    /\ providerRegistered' = FALSE
    /\ providerPending' = FALSE
    /\ pumpPending' =
        (demandPending \/ importanceCensusPending \/ qualityPending \/
         deadlineRecoveryPending \/
         (exactFrameRequired /\ ~exactFrameQueued))
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    demandPending,
                    demandCursorActive, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

CompleteProvider ==
    /\ providerPending
    /\ providerPending' = FALSE
    /\ demandPending' = automatic
    /\ demandCursorActive' = FALSE
    /\ qualityPending' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = FALSE
    /\ deadlineRecoveryPending' = FALSE
    /\ deadlineRetryAvailable' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    demandInterruptions, importanceCensusPending, renderState,
                    renderTargetRevision, renderTargetSemanticRevision,
                    retainedPresentation>>

CameraChange ==
    /\ viewRevision < MaxViewRevision
    /\ viewRevision' = viewRevision + 1
    /\ demandPending' = automatic
    /\ demandCursorActive' = FALSE
    /\ importanceCensusPending' = FALSE
    /\ qualityPending' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = FALSE
    /\ deadlineRecoveryPending' = FALSE
    /\ deadlineRetryAvailable' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandInterruptions,
                    renderState, renderTargetRevision,
                    renderTargetSemanticRevision,
                    retainedPresentation>>

TogglePolicy ==
    /\ policyChanges < MaxPolicyChanges
    /\ policyChanges' = policyChanges + 1
    /\ automatic' = ~automatic
    /\ demandPending' = ~automatic
    /\ demandCursorActive' = FALSE
    /\ importanceCensusPending' = FALSE
    /\ qualityPending' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = ~automatic
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = FALSE
    /\ deadlineRecoveryPending' = FALSE
    /\ deadlineRetryAvailable' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ UNCHANGED <<importanceRequests, semanticChanges,
                    viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered,
                    providersClosed, providerPending, renderState,
                    demandInterruptions,
                    renderTargetRevision, renderTargetSemanticRevision,
                    retainedPresentation>>

RequestImportanceCensus ==
    /\ automatic
    /\ importanceRequests < MaxImportanceRequests
    /\ ~providerPending
    /\ ~demandPending
    /\ ~qualityPending
    /\ ~importanceCensusPending
    /\ importanceRequests' = importanceRequests + 1
    /\ demandPending' = TRUE
    /\ demandCursorActive' = FALSE
    /\ importanceCensusPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = FALSE
    /\ frameSuperseded' = (renderState # "idle")
    /\ UNCHANGED <<automatic, policyChanges, semanticChanges,
                    viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, qualityPending, demandInterruptions,
                    renderState,
                    renderTargetRevision, renderTargetSemanticRevision,
                    deadlineRecoveryPending,
                    deadlineRetryAvailable, retainedPresentation>>

\* Selection, highlighting, and manipulator style changes are semantic-only.
\* They supersede an older framebuffer but preserve every LoD control fact and
\* capacity certificate.  A coalesced pre-existing capacity frame remains
\* capacity-relevant; this action never manufactures one.
SemanticMutation ==
    /\ semanticChanges < MaxSemanticChanges
    /\ semanticChanges' = semanticChanges + 1
    /\ semanticRevision' = semanticRevision + 1
    /\ exactFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    viewRevision, committedRevision,
                    committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, demandCursorActive, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    capacityFrameRequired, renderState,
                    renderTargetRevision, renderTargetSemanticRevision,
                    constraintWitness, deadlineRecoveryPending,
                    deadlineRetryAvailable, retainedPresentation>>

\* Demand is a level-triggered semantic obligation; the cursor is the finite
\* executable owner which can consume it.  A policy/source edge normally
\* starts this pass immediately in production, but the host must also be able
\* to recreate it later from the standing obligation alone.
StartDemand ==
    /\ automatic
    /\ demandPending
    /\ ~demandCursorActive
    /\ demandCursorActive' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

\* A stronger finite presentation owner may retire an in-progress ordinary
\* cursor without discharging its semantic demand.  Bound those interruptions
\* so the model checks the recovery seam without admitting an adversary which
\* can prevent every finite pass from completing forever.
InterruptDemand ==
    /\ demandCursorActive
    /\ demandInterruptions < MaxDemandInterruptions
    /\ demandCursorActive' = FALSE
    /\ demandInterruptions' = demandInterruptions + 1
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, importanceCensusPending, qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

CompleteDemand ==
    /\ automatic
    /\ demandPending
    /\ demandCursorActive
    /\ demandPending' = FALSE
    /\ demandCursorActive' = FALSE
    /\ importanceCensusPending' = FALSE
    /\ qualityPending' = TRUE
    /\ pumpPending' = TRUE
    /\ deadlineRecoveryPending' = FALSE
    /\ deadlineRetryAvailable' = TRUE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandInterruptions,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, renderState, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    retainedPresentation>>

CompleteQuality ==
    /\ automatic
    /\ qualityPending
    /\ qualityPending' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ deadlineRecoveryPending' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandPending, demandCursorActive,
                    demandInterruptions,
                    importanceCensusPending, renderState,
                    renderTargetRevision, renderTargetSemanticRevision,
                    constraintWitness,
                    deadlineRetryAvailable,
                    retainedPresentation>>

ConstrainQuality ==
    /\ automatic
    /\ qualityPending
    /\ qualityPending' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ capacityFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ constraintWitness' = TRUE
    /\ frameSuperseded' = (renderState # "idle")
    /\ deadlineRecoveryPending' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandPending, demandCursorActive,
                    demandInterruptions,
                    importanceCensusPending, renderState,
                    renderTargetRevision, renderTargetSemanticRevision,
                    deadlineRetryAvailable,
                    retainedPresentation>>

QueueExactFrame ==
    /\ exactFrameRequired
    /\ ~exactFrameQueued
    /\ ~providerPending
    /\ ~demandPending
    /\ ~qualityPending
    /\ renderState = "idle"
    /\ exactFrameQueued' = TRUE
    /\ renderState' = "pending"
    /\ renderTargetRevision' = viewRevision
    /\ renderTargetSemanticRevision' = semanticRevision
    /\ frameSuperseded' = FALSE
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandPending, demandCursorActive,
                    demandInterruptions, importanceCensusPending,
                    qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

BeginRender ==
    /\ renderState = "pending"
    /\ renderState' = "inflight"
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandPending, demandCursorActive,
                    demandInterruptions, importanceCensusPending,
                    qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    exactFrameQueued, pumpPending, renderTargetRevision,
                    renderTargetSemanticRevision, frameSuperseded,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

CompleteCurrentRender ==
    /\ renderState = "inflight"
    /\ renderTargetRevision = viewRevision
    /\ renderTargetSemanticRevision = semanticRevision
    /\ ~frameSuperseded
    /\ renderState' = "idle"
    /\ committedRevision' = viewRevision
    /\ committedSemanticRevision' = semanticRevision
    /\ exactFrameRequired' = FALSE
    /\ capacityFrameRequired' = FALSE
    /\ exactFrameQueued' = FALSE
    /\ frameSuperseded' = FALSE
    /\ pumpPending' =
        (providerPending \/ demandPending \/ importanceCensusPending \/
         qualityPending \/ deadlineRecoveryPending)
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, semanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, demandCursorActive, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    renderTargetRevision, renderTargetSemanticRevision,
                    constraintWitness, deadlineRecoveryPending,
                    deadlineRetryAvailable, retainedPresentation>>

CompleteCapacityDeadlineMiss ==
    /\ automatic
    /\ capacityFrameRequired
    /\ renderState = "inflight"
    /\ renderTargetRevision = viewRevision
    /\ renderTargetSemanticRevision = semanticRevision
    /\ ~frameSuperseded
    /\ deadlineRetryAvailable
    /\ renderState' = "idle"
    /\ qualityPending' = TRUE
    /\ exactFrameRequired' = TRUE
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ frameSuperseded' = FALSE
    /\ constraintWitness' = FALSE
    /\ deadlineRecoveryPending' = TRUE
    /\ deadlineRetryAvailable' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, demandCursorActive, demandInterruptions,
                    importanceCensusPending, capacityFrameRequired,
                    renderTargetRevision, renderTargetSemanticRevision,
                    retainedPresentation>>

\* Missing a semantic-only presentation deadline cannot manufacture a LoD
\* quality owner.  The exact style frame remains level-triggered and retries
\* once through the ordinary host render path.
CompleteSemanticDeadlineMiss ==
    /\ automatic
    /\ ~capacityFrameRequired
    /\ renderState = "inflight"
    /\ renderTargetRevision = viewRevision
    /\ renderTargetSemanticRevision = semanticRevision
    /\ ~frameSuperseded
    /\ deadlineRetryAvailable
    /\ renderState' = "idle"
    /\ exactFrameQueued' = FALSE
    /\ pumpPending' = TRUE
    /\ frameSuperseded' = FALSE
    /\ deadlineRetryAvailable' = FALSE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed, providerPending,
                    demandPending, demandCursorActive, demandInterruptions,
                    importanceCensusPending, qualityPending,
                    exactFrameRequired, capacityFrameRequired,
                    renderTargetRevision, renderTargetSemanticRevision,
                    constraintWitness, deadlineRecoveryPending,
                    retainedPresentation>>

CompleteStaleRender ==
    /\ renderState = "inflight"
    /\ (renderTargetRevision # viewRevision \/
        renderTargetSemanticRevision # semanticRevision \/ frameSuperseded)
    /\ renderState' = "idle"
    /\ exactFrameQueued' = FALSE
    /\ exactFrameRequired' = TRUE
    /\ frameSuperseded' = FALSE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<automatic, policyChanges, importanceRequests,
                    semanticChanges, viewRevision, committedRevision,
                    semanticRevision, committedSemanticRevision,
                    providerRegistered, providersClosed,
                    providerPending, demandPending, demandCursorActive,
                    demandInterruptions, importanceCensusPending,
                    qualityPending, capacityFrameRequired,
                    renderTargetRevision, renderTargetSemanticRevision,
                    constraintWitness,
                    deadlineRecoveryPending, deadlineRetryAvailable,
                    retainedPresentation>>

Next ==
    \/ RegisterProvider
    \/ PublishProviderWork
    \/ CloseProviders
    \/ CompleteProvider
    \/ CameraChange
    \/ TogglePolicy
    \/ RequestImportanceCensus
    \/ SemanticMutation
    \/ StartDemand
    \/ InterruptDemand
    \/ CompleteDemand
    \/ CompleteQuality
    \/ ConstrainQuality
    \/ QueueExactFrame
    \/ BeginRender
    \/ CompleteCapacityDeadlineMiss
    \/ CompleteSemanticDeadlineMiss
    \/ CompleteCurrentRender
    \/ CompleteStaleRender

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(CloseProviders)
    /\ WF_vars(CameraChange)
    /\ WF_vars(TogglePolicy)
    /\ WF_vars(RequestImportanceCensus)
    /\ WF_vars(SemanticMutation)
    /\ WF_vars(CompleteProvider)
    /\ WF_vars(StartDemand)
    /\ WF_vars(CompleteDemand)
    /\ WF_vars(CompleteQuality \/ ConstrainQuality)
    /\ WF_vars(QueueExactFrame)
    /\ WF_vars(BeginRender)
    /\ WF_vars(CompleteCapacityDeadlineMiss)
    /\ WF_vars(CompleteSemanticDeadlineMiss)
    /\ WF_vars(CompleteCurrentRender)
    /\ WF_vars(CompleteStaleRender)

PumpMatchesStandingWork == pumpPending = PumpReasons

DeadlineMissHasRecoveryOwner ==
    deadlineRecoveryPending => qualityPending /\ pumpPending

PresentationOnlyHasNoLodRecovery ==
    ~capacityFrameRequired =>
        ~qualityPending /\ ~deadlineRecoveryPending

DisabledOwnsNoAutomaticWork ==
    ~automatic =>
        ~demandPending /\ ~importanceCensusPending /\ ~qualityPending

RegisteredProviderIsNotWork ==
    providerRegistered /\ ~providerPending /\
    ~demandPending /\ ~importanceCensusPending /\ ~qualityPending /\
    ~(exactFrameRequired /\ ~exactFrameQueued)
        => ~pumpPending

ImportanceCensusHasDemandProducer ==
    importanceCensusPending =>
        demandPending /\ (demandCursorActive \/ pumpPending)

DemandCursorIsOwned == demandCursorActive => automatic /\ demandPending

TerminalIsCurrentAndQuiescent ==
    Outcome \in {"ready", "constrained"} =>
        /\ committedRevision = viewRevision
        /\ committedSemanticRevision = semanticRevision
        /\ ~ForegroundWork
        /\ ~exactFrameRequired
        /\ ~capacityFrameRequired

ConstrainedHasWitness == Outcome = "constrained" => constraintWitness

ReadyHasNoConstraint == Outcome = "ready" => ~constraintWitness

RetainedPresentationSurvives == retainedPresentation

EventuallySettlesAfterInputsClose ==
    []((viewRevision = MaxViewRevision /\
        policyChanges = MaxPolicyChanges /\
        importanceRequests = MaxImportanceRequests /\
        semanticChanges = MaxSemanticChanges /\ providersClosed) =>
       <> (Outcome \in {"ready", "constrained"}))

=============================================================================
