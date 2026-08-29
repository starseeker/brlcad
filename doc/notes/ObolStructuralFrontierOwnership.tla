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

CONSTANT MaxPointSteps

Frontiers == {"clear", "aggregatable", "meshRequired"}
Owners == {"none", "admissionFrame", "structuralRepair",
           "structuralPrefetch", "pointPublication", "stableCalibration"}

VARIABLES frontier, owner, pointStep, pointSettled, structuralBacking,
          candidatePending, candidateStep, capacityRevision, terminal

vars == <<frontier, owner, pointStep, pointSettled, structuralBacking,
          candidatePending, candidateStep, capacityRevision, terminal>>

TypeOK ==
    /\ frontier \in Frontiers
    /\ owner \in Owners
    /\ pointStep \in 0..MaxPointSteps
    /\ pointSettled \in BOOLEAN
    /\ structuralBacking \in BOOLEAN
    /\ candidatePending \in BOOLEAN
    /\ candidateStep \in 0..MaxPointSteps
    /\ capacityRevision \in 0..MaxPointSteps
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
    /\ terminal = FALSE

ClaimAdmissionFrame ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "aggregatable"
    /\ pointStep < MaxPointSteps
    /\ owner' = "admissionFrame"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    terminal>>

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
                    terminal>>

ClaimStructuralRepair ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "meshRequired"
    /\ owner' = "structuralRepair"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    terminal>>

CompleteStructuralRepair ==
    /\ ~terminal
    /\ owner = "structuralRepair"
    /\ frontier' = "clear"
    /\ owner' = "none"
    /\ pointSettled' = (pointStep = 0)
    /\ structuralBacking' \in IF pointStep = 0 THEN {FALSE}
                             ELSE BOOLEAN
    /\ UNCHANGED <<pointStep, candidatePending, candidateStep,
                    capacityRevision, terminal>>

ClaimStableCalibration ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "clear"
    /\ pointStep > 0
    /\ ~pointSettled
    /\ owner' = "stableCalibration"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    candidatePending, candidateStep, capacityRevision,
                    terminal>>

CompleteStableCalibration ==
    /\ ~terminal
    /\ owner = "stableCalibration"
    /\ \/ /\ owner' = "none"
           /\ pointSettled' = TRUE
           /\ UNCHANGED <<frontier, pointStep, structuralBacking,
                           candidatePending, candidateStep,
                           capacityRevision, terminal>>
       \/ /\ pointStep > 0
           /\ ~structuralBacking
           /\ owner' = "none"
           /\ pointStep' = pointStep - 1
           /\ pointSettled' = (pointStep' = 0)
           /\ capacityRevision' = capacityRevision + 1
           /\ UNCHANGED <<frontier, structuralBacking, candidatePending,
                           candidateStep, terminal>>
       \/ /\ pointStep > 0
           /\ structuralBacking
           /\ owner' = "structuralPrefetch"
           /\ candidatePending' = TRUE
           /\ candidateStep' = pointStep - 1
           /\ UNCHANGED <<frontier, pointStep, pointSettled,
                           structuralBacking, capacityRevision, terminal>>

CompleteStructuralPrefetch ==
    /\ ~terminal
    /\ owner = "structuralPrefetch"
    /\ pointStep > 0
    /\ candidatePending
    /\ candidateStep = pointStep - 1
    /\ \/ /\ owner' = "pointPublication"
           /\ UNCHANGED <<pointSettled, structuralBacking,
                           candidatePending, candidateStep>>
       \/ /\ owner' = "none"
           /\ pointSettled' = TRUE
           /\ candidatePending' = FALSE
           /\ candidateStep' = 0
           /\ UNCHANGED structuralBacking
    /\ UNCHANGED <<frontier, pointStep, capacityRevision, terminal>>

\* The old point cut remains active throughout preload.  Only this exact
\* publication edge may expose the candidate and invalidate the capacity
\* certificate.  If the exact frame still cannot replace its structural
\* frontier, rejection preserves both the old cut and revision.
CompletePointPublication ==
    /\ ~terminal
    /\ owner = "pointPublication"
    /\ candidatePending
    /\ candidateStep = pointStep - 1
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
    /\ UNCHANGED <<frontier, terminal>>

PublishReady ==
    /\ ~terminal
    /\ frontier = "clear"
    /\ owner = "none"
    /\ ~candidatePending
    /\ (pointStep = 0 \/ pointSettled)
    /\ terminal' = TRUE
    /\ UNCHANGED <<frontier, owner, pointStep, pointSettled,
                    structuralBacking, candidatePending, candidateStep,
                    capacityRevision>>

ClaimSuccessor == ClaimAdmissionFrame \/ ClaimStructuralRepair \/
                  ClaimStableCalibration
CompleteSuccessor == CompleteAdmissionFrame \/ CompleteStructuralRepair \/
                     CompleteStableCalibration \/
                     CompleteStructuralPrefetch \/
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
    /\ (pointStep = 0 \/ pointSettled)

OwnerMatchesFrontier ==
    /\ (owner = "admissionFrame" => frontier = "aggregatable")
    /\ (owner = "structuralRepair" => frontier = "meshRequired")
    /\ (owner = "structuralPrefetch" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled /\
           structuralBacking /\ candidatePending /\
           candidateStep = pointStep - 1)
    /\ (owner = "pointPublication" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled /\
           candidatePending /\ candidateStep = pointStep - 1)
    /\ (owner = "stableCalibration" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled)

PrivateCandidateHasOneOwner ==
    candidatePending =>
        /\ owner \in {"structuralPrefetch", "pointPublication"}
        /\ candidateStep = pointStep - 1

\* Preload completion never mutates the public cut or capacity revision.
\* The sole candidate-publication transition either advances both atomically
\* or clears the candidate while preserving both.
PrivatePreloadCannotPublish ==
    owner = "structuralPrefetch" =>
        /\ pointStep > 0
        /\ candidateStep = pointStep - 1

OwnerlessFrontierHasSuccessor ==
    owner = "none" /\
    (frontier # "clear" \/ (pointStep > 0 /\ ~pointSettled)) =>
        ENABLED ClaimSuccessor

EventuallyReady == <>terminal

=============================================================================
