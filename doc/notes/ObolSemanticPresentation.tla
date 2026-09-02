-------------------- MODULE ObolSemanticPresentation --------------------
\* Exact-frame ownership for selection, highlighting, and other retained CAD
\* style mutations.  These changes do not alter geometry allocation or its
\* capacity certificate or advance LoD planning, but readiness must wait until
\* a frame which began after the newest semantic revision has completed.

EXTENDS Naturals, TLC

CONSTANT MaxSemanticRevision

ASSUME MaxSemanticRevision > 0

FrameStates == {"idle", "pending", "inflight"}

VARIABLES semanticRevision, committedRevision, frameRevision, frameState,
          requestRequired, recoveryFaultAvailable, frameCapacityRelevant,
          framePlanningRelevant

vars == <<semanticRevision, committedRevision, frameRevision, frameState,
          requestRequired, recoveryFaultAvailable, frameCapacityRelevant,
          framePlanningRelevant>>

Init ==
    /\ semanticRevision = 0
    /\ committedRevision = 0
    /\ frameRevision = 0
    /\ frameState = "idle"
    /\ requestRequired = FALSE
    /\ recoveryFaultAvailable = TRUE
    /\ frameCapacityRelevant = FALSE
    /\ framePlanningRelevant = FALSE

MutateStyle ==
    /\ semanticRevision < MaxSemanticRevision
    /\ semanticRevision' = semanticRevision + 1
    /\ requestRequired' = TRUE
    /\ UNCHANGED <<committedRevision, frameRevision, frameState,
                    recoveryFaultAvailable, frameCapacityRelevant,
                    framePlanningRelevant>>

\* Abstract the implementation seam which motivates idle recovery: retained
\* CAD state can become stale without its original host request remaining.
\* The fault is bounded so liveness can require the recovery transition.
LoseOwner ==
    /\ recoveryFaultAvailable
    /\ requestRequired
    /\ frameState = "idle"
    /\ requestRequired' = FALSE
    /\ recoveryFaultAvailable' = FALSE
    /\ UNCHANGED <<semanticRevision, committedRevision, frameRevision,
                    frameState, frameCapacityRelevant,
                    framePlanningRelevant>>

OwnerlessStale ==
    /\ committedRevision # semanticRevision
    /\ frameState = "idle"
    /\ ~requestRequired

RecoverOwnerless ==
    /\ OwnerlessStale
    /\ requestRequired' = TRUE
    /\ UNCHANGED <<semanticRevision, committedRevision, frameRevision,
                    frameState, recoveryFaultAvailable,
                    frameCapacityRelevant, framePlanningRelevant>>

QueueFrame ==
    /\ requestRequired
    /\ frameState = "idle"
    /\ frameRevision' = semanticRevision
    /\ frameState' = "pending"
    /\ requestRequired' = FALSE
    /\ frameCapacityRelevant' = FALSE
    /\ framePlanningRelevant' = FALSE
    /\ UNCHANGED <<semanticRevision, committedRevision,
                    recoveryFaultAvailable>>

BeginFrame ==
    /\ frameState = "pending"
    /\ frameState' = "inflight"
    /\ UNCHANGED <<semanticRevision, committedRevision, frameRevision,
                    requestRequired, recoveryFaultAvailable,
                    frameCapacityRelevant, framePlanningRelevant>>

CompleteCurrentFrame ==
    /\ frameState = "inflight"
    /\ frameRevision = semanticRevision
    /\ committedRevision' = semanticRevision
    /\ frameState' = "idle"
    /\ frameCapacityRelevant' = FALSE
    /\ framePlanningRelevant' = FALSE
    /\ UNCHANGED <<semanticRevision, frameRevision, requestRequired,
                    recoveryFaultAvailable>>

CompleteStaleFrame ==
    /\ frameState = "inflight"
    /\ frameRevision # semanticRevision
    /\ frameState' = "idle"
    /\ requestRequired' = TRUE
    /\ frameCapacityRelevant' = FALSE
    /\ framePlanningRelevant' = FALSE
    /\ UNCHANGED <<semanticRevision, committedRevision, frameRevision,
                    recoveryFaultAvailable>>

Next == MutateStyle \/ LoseOwner \/ RecoverOwnerless \/ QueueFrame \/ BeginFrame \/
        CompleteCurrentFrame \/ CompleteStaleFrame

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(QueueFrame)
    /\ WF_vars(RecoverOwnerless)
    /\ WF_vars(BeginFrame)
    /\ WF_vars(CompleteCurrentFrame)
    /\ WF_vars(CompleteStaleFrame)

TypeOK ==
    /\ semanticRevision \in 0..MaxSemanticRevision
    /\ committedRevision \in 0..MaxSemanticRevision
    /\ committedRevision <= semanticRevision
    /\ frameRevision \in 0..MaxSemanticRevision
    /\ frameState \in FrameStates
    /\ requestRequired \in BOOLEAN
    /\ recoveryFaultAvailable \in BOOLEAN
    /\ frameCapacityRelevant \in BOOLEAN
    /\ framePlanningRelevant \in BOOLEAN

Ready == committedRevision = semanticRevision /\ frameState = "idle" /\
         ~requestRequired

StaleFrameCannotCommit == committedRevision <= frameRevision

\* Selection/style repair is a presentation transaction.  A pre-existing LoD
\* capacity frame may absorb it in the implementation, but this isolated
\* semantic owner must never manufacture capacity or planning evidence.
SemanticFrameIsPresentationOnly ==
    frameState = "idle" \/
        (~frameCapacityRelevant /\ ~framePlanningRelevant)

OwnedSuccessorNeedsNoRecovery ==
    (requestRequired \/ frameState # "idle") => ~OwnerlessStale

EventuallyReadyAfterMutationsStop ==
    [](semanticRevision = MaxSemanticRevision => <>Ready)

=============================================================================
