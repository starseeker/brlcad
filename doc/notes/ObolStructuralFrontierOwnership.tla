-------------------- MODULE ObolStructuralFrontierOwnership --------------------
\* A focused model of the quiet structural-presentation frontier.  An exact
\* framebuffer may prove that structural fallbacks remain after every other
\* producer has gone idle.  That state must own one of two finite successors:
\* another point-aggregation frame, or mesh repair after the point policy is
\* exhausted.  It may never publish the terminal "View ready" promise.

EXTENDS Naturals

CONSTANT MaxPointSteps

Frontiers == {"clear", "aggregatable", "meshRequired"}
Owners == {"none", "pointFrame", "structuralRepair"}

VARIABLES frontier, owner, pointStep, terminal

vars == <<frontier, owner, pointStep, terminal>>

TypeOK ==
    /\ frontier \in Frontiers
    /\ owner \in Owners
    /\ pointStep \in 0..MaxPointSteps
    /\ terminal \in BOOLEAN

Init ==
    /\ frontier \in Frontiers
    /\ owner = "none"
    /\ pointStep = 0
    /\ terminal = FALSE

ClaimPointFrame ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "aggregatable"
    /\ pointStep < MaxPointSteps
    /\ owner' = "pointFrame"
    /\ UNCHANGED <<frontier, pointStep, terminal>>

CompletePointFrame ==
    /\ ~terminal
    /\ owner = "pointFrame"
    /\ owner' = "none"
    /\ pointStep' = pointStep + 1
    /\ frontier' \in
        IF pointStep' = MaxPointSteps
        THEN {"clear", "meshRequired"}
        ELSE {"clear", "aggregatable"}
    /\ UNCHANGED terminal

ClaimStructuralRepair ==
    /\ ~terminal
    /\ owner = "none"
    /\ frontier = "meshRequired"
    /\ owner' = "structuralRepair"
    /\ UNCHANGED <<frontier, pointStep, terminal>>

CompleteStructuralRepair ==
    /\ ~terminal
    /\ owner = "structuralRepair"
    /\ frontier' = "clear"
    /\ owner' = "none"
    /\ UNCHANGED <<pointStep, terminal>>

PublishReady ==
    /\ ~terminal
    /\ frontier = "clear"
    /\ owner = "none"
    /\ terminal' = TRUE
    /\ UNCHANGED <<frontier, owner, pointStep>>

ClaimSuccessor == ClaimPointFrame \/ ClaimStructuralRepair
CompleteSuccessor == CompletePointFrame \/ CompleteStructuralRepair

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

NoPrematureReady == terminal => frontier = "clear" /\ owner = "none"

OwnerMatchesFrontier ==
    /\ (owner = "pointFrame" => frontier = "aggregatable")
    /\ (owner = "structuralRepair" => frontier = "meshRequired")

OwnerlessFrontierHasSuccessor ==
    frontier # "clear" /\ owner = "none" => ENABLED ClaimSuccessor

EventuallyReady == <>terminal

=============================================================================
