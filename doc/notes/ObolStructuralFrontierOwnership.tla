-------------------- MODULE ObolStructuralFrontierOwnership --------------------
\* A focused model of the quiet structural-presentation frontier.  An exact
\* framebuffer may prove that structural fallbacks remain after every other
\* producer has gone idle.  That state must own one of two finite successors:
\* another structural-admission classifier frame, or mesh repair after the
\* point policy is exhausted.  Once repair clears the frontier, a retained
\* coarse point cut also requires one stable-quality frame before it may
\* publish the terminal "View ready" promise.  A finer point cut which would
\* expose hidden structural records first preloads that exact sparse frontier;
\* only its completed publication may commit the finer cut.  Structural
\* admission, sparse prefetch, and stable calibration are distinct owners.

EXTENDS Naturals

CONSTANT MaxPointSteps

Frontiers == {"clear", "aggregatable", "meshRequired"}
Owners == {"none", "admissionFrame", "structuralRepair",
           "structuralPrefetch", "stableCalibration"}

VARIABLES frontier, owner, pointStep, pointSettled, structuralBacking,
          terminal

vars == <<frontier, owner, pointStep, pointSettled, structuralBacking,
          terminal>>

TypeOK ==
    /\ frontier \in Frontiers
    /\ owner \in Owners
    /\ pointStep \in 0..MaxPointSteps
    /\ pointSettled \in BOOLEAN
    /\ structuralBacking \in BOOLEAN
    /\ terminal \in BOOLEAN

Init ==
    /\ frontier \in Frontiers
    /\ owner = "none"
    /\ pointStep = 0
    /\ pointSettled = TRUE
    /\ structuralBacking = FALSE
    /\ terminal = FALSE

ClaimAdmissionFrame ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "aggregatable"
    /\ pointStep < MaxPointSteps
    /\ owner' = "admissionFrame"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
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
    /\ UNCHANGED terminal

ClaimStructuralRepair ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "meshRequired"
    /\ owner' = "structuralRepair"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    terminal>>

CompleteStructuralRepair ==
    /\ ~terminal
    /\ owner = "structuralRepair"
    /\ frontier' = "clear"
    /\ owner' = "none"
    /\ pointSettled' = (pointStep = 0)
    /\ structuralBacking' \in IF pointStep = 0 THEN {FALSE}
                             ELSE BOOLEAN
    /\ UNCHANGED <<pointStep, terminal>>

ClaimStableCalibration ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "clear"
    /\ pointStep > 0
    /\ ~pointSettled
    /\ owner' = "stableCalibration"
    /\ UNCHANGED <<frontier, pointStep, pointSettled, structuralBacking,
                    terminal>>

CompleteStableCalibration ==
    /\ ~terminal
    /\ owner = "stableCalibration"
    /\ \/ /\ owner' = "none"
           /\ pointSettled' = TRUE
           /\ UNCHANGED <<frontier, pointStep, structuralBacking,
                           terminal>>
       \/ /\ pointStep > 0
           /\ ~structuralBacking
           /\ owner' = "none"
           /\ pointStep' = pointStep - 1
           /\ pointSettled' = (pointStep' = 0)
           /\ UNCHANGED <<frontier, structuralBacking, terminal>>
       \/ /\ pointStep > 0
           /\ structuralBacking
           /\ owner' = "structuralPrefetch"
           /\ UNCHANGED <<frontier, pointStep, pointSettled,
                           structuralBacking, terminal>>

CompleteStructuralPrefetch ==
    /\ ~terminal
    /\ owner = "structuralPrefetch"
    /\ pointStep > 0
    /\ owner' = "none"
    /\ \/ /\ pointStep' = pointStep - 1
           /\ pointSettled' = (pointStep' = 0)
           /\ structuralBacking' \in BOOLEAN
       \/ /\ UNCHANGED <<pointStep, structuralBacking>>
           /\ pointSettled' = TRUE
    /\ UNCHANGED <<frontier, terminal>>

PublishReady ==
    /\ ~terminal
    /\ frontier = "clear"
    /\ owner = "none"
    /\ (pointStep = 0 \/ pointSettled)
    /\ terminal' = TRUE
    /\ UNCHANGED <<frontier, owner, pointStep, pointSettled,
                    structuralBacking>>

ClaimSuccessor == ClaimAdmissionFrame \/ ClaimStructuralRepair \/
                  ClaimStableCalibration
CompleteSuccessor == CompleteAdmissionFrame \/ CompleteStructuralRepair \/
                     CompleteStableCalibration \/
                     CompleteStructuralPrefetch

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
    /\ (pointStep = 0 \/ pointSettled)

OwnerMatchesFrontier ==
    /\ (owner = "admissionFrame" => frontier = "aggregatable")
    /\ (owner = "structuralRepair" => frontier = "meshRequired")
    /\ (owner = "structuralPrefetch" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled /\
           structuralBacking)
    /\ (owner = "stableCalibration" =>
           frontier = "clear" /\ pointStep > 0 /\ ~pointSettled)

OwnerlessFrontierHasSuccessor ==
    owner = "none" /\
    (frontier # "clear" \/ (pointStep > 0 /\ ~pointSettled)) =>
        ENABLED ClaimSuccessor

EventuallyReady == <>terminal

=============================================================================
