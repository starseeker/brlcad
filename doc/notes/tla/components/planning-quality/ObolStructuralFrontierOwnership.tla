-------------------- MODULE ObolStructuralFrontierOwnership --------------------
\* A focused model of the quiet structural-presentation frontier.  An exact
\* framebuffer may prove that structural fallbacks remain after every other
\* producer has gone idle.  That state must own one of two finite successors:
\* another structural-admission classifier frame, or mesh repair after the
\* point policy is exhausted.  Once repair clears the frontier, a retained
\* coarse point cut also requires one stable-quality frame before it may
\* publish the terminal "View ready" promise.  A finer point cut which would
\* expose hidden structural records first preloads that exact sparse frontier
\* while the old cut and capacity certificate remain public.  The candidate
\* then commits exactly once through an exact publication, or rejects exactly
\* once without changing either public fact.  Structural admission, sparse
\* prefetch, publication, and stable calibration are distinct owners.

EXTENDS Naturals

CONSTANTS MaxPointSteps, MaxPrefetch, MaxPrefetchBatch

Frontiers == {"clear", "aggregatable", "meshRequired"}
Owners == {"none", "admissionFrame", "structuralRepair",
           "terminalProxy", "structuralPrefetch", "pointPublication",
           "stableCalibration"}

VARIABLES frontier, owner, pointStep, pointSettled, structuralBacking,
          candidatePending, candidateStep, capacityRevision,
          structuralCapacityRejected, prefetchRemaining, prefetchBatch,
          terminal

vars == <<frontier, owner, pointStep, pointSettled, structuralBacking,
          candidatePending, candidateStep, capacityRevision,
          structuralCapacityRejected, prefetchRemaining, prefetchBatch,
          terminal>>

TypeOK ==
    /\ frontier \in Frontiers
    /\ owner \in Owners
    /\ pointStep \in 0..MaxPointSteps
    /\ pointSettled \in BOOLEAN
    /\ structuralBacking \in BOOLEAN
    /\ candidatePending \in BOOLEAN
    /\ candidateStep \in 0..MaxPointSteps
    /\ capacityRevision \in 0..MaxPointSteps
    /\ structuralCapacityRejected \in BOOLEAN
    /\ prefetchRemaining \in 0..MaxPrefetch
    /\ prefetchBatch \in 0..MaxPrefetchBatch
    /\ terminal \in BOOLEAN

Init ==
    /\ frontier \in Frontiers
    /\ owner = "none"
    /\ pointStep = 0
    /\ pointSettled = TRUE
    /\ structuralBacking = FALSE
    /\ candidatePending = FALSE
    /\ candidateStep = 0
    /\ capacityRevision = 0
    /\ structuralCapacityRejected = FALSE
    /\ prefetchRemaining = 0
    /\ prefetchBatch = 0
    /\ terminal = FALSE

ClaimAdmissionFrame ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "aggregatable"
    /\ pointStep < MaxPointSteps
    /\ owner' = "admissionFrame"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, prefetchRemaining,
                    prefetchBatch, terminal>>

CompleteAdmissionFrame ==
    /\ ~terminal
    /\ owner = "admissionFrame"
    /\ owner' = "none"
    /\ pointStep' = pointStep + 1
    /\ frontier' \in
        IF pointStep' = MaxPointSteps
        THEN {"clear", "meshRequired"}
        ELSE {"clear", "aggregatable"}
    /\ pointSettled' = FALSE
    /\ structuralBacking' = TRUE
    /\ UNCHANGED <<candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, prefetchRemaining,
                    prefetchBatch, terminal>>

ClaimStructuralRepair ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "meshRequired"
    /\ ~structuralCapacityRejected
    /\ owner' = "structuralRepair"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, prefetchRemaining,
                    prefetchBatch, terminal>>

CompleteStructuralRepairSuccess ==
    /\ ~terminal
    /\ owner = "structuralRepair"
    /\ owner' = "none"
    /\ frontier' = "clear"
    /\ pointSettled' = (pointStep = 0)
    /\ structuralBacking' \in IF pointStep = 0 THEN {FALSE}
                             ELSE BOOLEAN
    /\ structuralCapacityRejected' = structuralCapacityRejected
    /\ UNCHANGED <<pointStep, candidatePending, candidateStep,
                    capacityRevision, prefetchRemaining, prefetchBatch,
                    terminal>>

RejectStructuralRepair ==
    /\ ~terminal
    /\ owner = "structuralRepair"
    /\ owner' = "none"
    /\ structuralCapacityRejected' = TRUE
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    prefetchRemaining, prefetchBatch, terminal>>

ClaimTerminalProxy ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "meshRequired"
    /\ structuralCapacityRejected
    /\ owner' = "terminalProxy"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, prefetchRemaining,
                    prefetchBatch, terminal>>

CompleteTerminalProxy ==
    /\ ~terminal
    /\ owner = "terminalProxy"
    /\ frontier' = "clear"
    /\ owner' = "none"
    /\ pointSettled' = TRUE
    /\ structuralBacking' = FALSE
    /\ UNCHANGED <<pointStep, candidatePending, candidateStep,
                    capacityRevision, structuralCapacityRejected,
                    prefetchRemaining, prefetchBatch, terminal>>

ClaimStableCalibration ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "clear"
    /\ pointStep > 0
    /\ ~pointSettled
    /\ owner' = "stableCalibration"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, prefetchRemaining,
                    prefetchBatch, terminal>>

CompleteStableCalibration ==
    /\ ~terminal
    /\ owner = "stableCalibration"
    /\ \/ /\ owner' = "none"
           /\ pointSettled' = TRUE
           /\ UNCHANGED <<frontier, pointStep, structuralBacking,
                           candidatePending, candidateStep,
                           capacityRevision, structuralCapacityRejected,
                           prefetchRemaining, prefetchBatch, terminal>>
       \/ /\ pointStep > 0
           /\ ~structuralBacking
           /\ owner' = "none"
           /\ pointStep' = pointStep - 1
           /\ pointSettled' = (pointStep' = 0)
           /\ capacityRevision' = capacityRevision + 1
           /\ UNCHANGED <<frontier, structuralBacking, candidatePending,
                           candidateStep, structuralCapacityRejected,
                           prefetchRemaining, prefetchBatch, terminal>>
       \/ /\ pointStep > 0
           /\ structuralBacking
           /\ owner' = "structuralPrefetch"
           /\ candidatePending' = TRUE
           /\ candidateStep' = pointStep - 1
           /\ prefetchRemaining' \in 1..MaxPrefetch
           /\ prefetchBatch' \in 1..MaxPrefetchBatch
           /\ prefetchBatch' <= prefetchRemaining'
           /\ UNCHANGED <<frontier, pointStep, pointSettled,
                           structuralBacking, capacityRevision,
                           structuralCapacityRejected, terminal>>

\* Each exact bounded-batch completion either clears the frontier or exposes a
\* strictly smaller remainder.  The completed batch can replace between one
\* and every selected occurrence; a nonzero remainder receives another bounded
\* batch.  Provider publication changes inventory while fulfilling this action,
\* but does not change the candidate's semantic occurrence/projection domain
\* and therefore cannot cancel it.
AdvanceStructuralPrefetch ==
    /\ ~terminal
    /\ owner = "structuralPrefetch"
    /\ pointStep > 0
    /\ candidatePending
    /\ candidateStep = pointStep - 1
    /\ prefetchRemaining > 0
    /\ prefetchBatch > 0
    /\ \E nextRemaining \in 0..(prefetchRemaining - 1):
          /\ prefetchRemaining <= nextRemaining + prefetchBatch
          /\ prefetchRemaining' = nextRemaining
          /\ IF nextRemaining = 0
                THEN /\ owner' = "pointPublication"
                     /\ prefetchBatch' = 0
                ELSE /\ owner' = "structuralPrefetch"
                     /\ prefetchBatch' \in 1..MaxPrefetchBatch
                     /\ prefetchBatch' <= nextRemaining
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    structuralCapacityRejected, terminal>>

\* A true domain change or a pass with no strict progress rejects the private
\* candidate without changing the public point cut or capacity evidence.
RejectStructuralPrefetch ==
    /\ ~terminal
    /\ owner = "structuralPrefetch"
    /\ candidatePending
    /\ prefetchRemaining > 0
    /\ owner' = "none"
    /\ pointSettled' = TRUE
    /\ candidatePending' = FALSE
    /\ candidateStep' = 0
    /\ prefetchRemaining' = 0
    /\ prefetchBatch' = 0
    /\ UNCHANGED <<frontier, pointStep, capacityRevision,
                    structuralBacking, structuralCapacityRejected, terminal>>

\* The old point cut remains active throughout preload.  Only this exact
\* publication edge may expose the candidate and invalidate the capacity
\* certificate.  If the exact frame still cannot replace its structural
\* frontier, rejection preserves both the old cut and revision.
CompletePointPublication ==
    /\ ~terminal
    /\ owner = "pointPublication"
    /\ candidatePending
    /\ candidateStep = pointStep - 1
    /\ prefetchRemaining = 0
    /\ prefetchBatch = 0
    /\ \/ /\ owner' = "none"
           /\ pointStep' = candidateStep
           /\ pointSettled' = (candidateStep = 0)
           /\ candidatePending' = FALSE
           /\ candidateStep' = 0
           /\ capacityRevision' = capacityRevision + 1
           /\ structuralBacking' \in BOOLEAN
       \/ /\ owner' = "none"
           /\ pointSettled' = TRUE
           /\ candidatePending' = FALSE
           /\ candidateStep' = 0
           /\ UNCHANGED <<pointStep, structuralBacking,
                           capacityRevision>>
    /\ UNCHANGED <<frontier, structuralCapacityRejected,
                    prefetchRemaining, prefetchBatch, terminal>>

PublishReady ==
    /\ ~terminal
    /\ frontier = "clear"
    /\ owner = "none"
    /\ ~candidatePending
    /\ prefetchRemaining = 0
    /\ prefetchBatch = 0
    /\ (pointStep = 0 \/ pointSettled)
    /\ terminal' = TRUE
    /\ UNCHANGED <<frontier, owner, pointStep, pointSettled,
                    structuralBacking, candidatePending, candidateStep,
                    capacityRevision, structuralCapacityRejected,
                    prefetchRemaining, prefetchBatch>>

ClaimSuccessor == ClaimAdmissionFrame \/ ClaimStructuralRepair \/
                  ClaimTerminalProxy \/ ClaimStableCalibration
CompleteSuccessor == CompleteAdmissionFrame \/
                     CompleteStructuralRepairSuccess \/
                     RejectStructuralRepair \/
                     CompleteTerminalProxy \/
                     CompleteStableCalibration \/
                     AdvanceStructuralPrefetch \/
                     RejectStructuralPrefetch \/
                     CompletePointPublication

Next ==
    \/ ClaimSuccessor
    \/ CompleteSuccessor
    \/ PublishReady

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ClaimSuccessor)
    /\ WF_vars(CompleteSuccessor)
    /\ WF_vars(PublishReady)

NoPrematureReady == terminal =>
    /\ frontier = "clear"
    /\ owner = "none"
    /\ ~candidatePending
    /\ prefetchRemaining = 0
    /\ prefetchBatch = 0
    /\ (pointStep = 0 \/ pointSettled)

OwnerMatchesFrontier ==
    /\ (owner = "admissionFrame" => frontier = "aggregatable")
    /\ (owner = "structuralRepair" => frontier = "meshRequired")
    /\ (owner = "terminalProxy" =>
           frontier = "meshRequired" /\ structuralCapacityRejected)
    /\ (owner = "structuralPrefetch" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled /\
           structuralBacking /\ candidatePending /\
           candidateStep = pointStep - 1 /\ prefetchRemaining > 0 /\
           prefetchBatch > 0 /\ prefetchBatch <= prefetchRemaining)
    /\ (owner = "pointPublication" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled /\
           candidatePending /\ candidateStep = pointStep - 1 /\
           prefetchRemaining = 0)
    /\ (owner = "stableCalibration" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled)

PrivateCandidateHasOneOwner ==
    candidatePending =>
        /\ owner \in {"structuralPrefetch", "pointPublication"}
        /\ candidateStep = pointStep - 1

PrefetchRankOwned ==
    /\ (prefetchRemaining > 0) <=> (owner = "structuralPrefetch")
    /\ (prefetchBatch > 0) <=> (owner = "structuralPrefetch")
    /\ prefetchBatch <= prefetchRemaining

BoundedPrefetchBatch == prefetchBatch <= MaxPrefetchBatch

\* Preload completion never mutates the public cut or capacity revision.
\* The sole candidate-publication transition either advances both atomically
\* or clears the candidate while preserving both.
PrivatePreloadCannotPublish ==
    owner = "structuralPrefetch" =>
        /\ pointStep > 0
        /\ candidateStep = pointStep - 1
        /\ prefetchRemaining > 0
        /\ prefetchBatch > 0

TerminalProxyRequiresCapacityWitness ==
    owner = "terminalProxy" => structuralCapacityRejected

OwnerlessFrontierHasSuccessor ==
    owner = "none" /\
    (frontier # "clear" \/ (pointStep > 0 /\ ~pointSettled)) =>
        ENABLED ClaimSuccessor

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => terminal

EventuallyReady == <>terminal

=============================================================================
