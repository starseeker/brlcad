----------------------- MODULE ObolSubmissionPass -----------------------
\* One bounded model of the submission cursor and its full-rescan debt.
\* Discovery may finish an exact append pass before compact inventory is
\* complete.  That pauses the cursor without consuming the rescan.  Fair
\* inventory publication and cursor ownership must eventually consume it.

EXTENDS TLC

States == {"idle", "active", "idle-rescan", "active-rescan"}

VARIABLES state, inventoryComplete, coveragePending, pumpPending

vars == <<state, inventoryComplete, coveragePending, pumpPending>>

Active == state \in {"active", "active-rescan"}
RescanPending == state \in {"idle-rescan", "active-rescan"}

Init ==
    /\ state \in {"idle", "idle-rescan"}
    /\ inventoryComplete = (state = "idle")
    /\ coveragePending = TRUE
    /\ pumpPending = (state = "idle")

PublishCompleteInventory ==
    /\ ~inventoryComplete
    /\ inventoryComplete' = TRUE
    /\ UNCHANGED <<state, coveragePending, pumpPending>>

ResumeRescan ==
    /\ state = "idle-rescan"
    /\ state' = "active-rescan"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending>>

FinishRescanPass ==
    /\ state = "active-rescan"
    /\ state' = IF inventoryComplete THEN "idle" ELSE "idle-rescan"
    /\ coveragePending' = IF inventoryComplete THEN FALSE
                           ELSE coveragePending
    /\ pumpPending' = FALSE
    /\ UNCHANGED inventoryComplete

\* A pose-continuity or stronger-owner pause may leave the semantic coverage
\* obligation active after its cursor retires.  Once producer inventory is
\* closed, the level-triggered owner must restore one bounded cursor without
\* waiting for another camera or source event.
ResumeCoverage ==
    /\ state = "idle"
    /\ inventoryComplete
    /\ coveragePending
    /\ pumpPending
    /\ state' = "active"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending>>

FinishCoveragePass ==
    /\ state = "active"
    /\ coveragePending
    /\ state' = "idle"
    /\ coveragePending' = FALSE
    /\ pumpPending' = FALSE
    /\ UNCHANGED inventoryComplete

Next == PublishCompleteInventory \/ ResumeRescan \/ FinishRescanPass \/
        ResumeCoverage \/ FinishCoveragePass

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishCompleteInventory)
    /\ WF_vars(ResumeRescan)
    /\ WF_vars(FinishRescanPass)
    /\ WF_vars(ResumeCoverage)
    /\ WF_vars(FinishCoveragePass)

TypeOK ==
    /\ state \in States
    /\ inventoryComplete \in BOOLEAN
    /\ coveragePending \in BOOLEAN
    /\ pumpPending \in BOOLEAN

PausedPassPreservesDebt == state = "idle-rescan" => ~Active /\ RescanPending

ActiveRescanOwnsBothFacts ==
    state = "active-rescan" => Active /\ RescanPending

EventuallyConsumesRescan == <> (state = "idle" /\ ~RescanPending)

CoverageHasProducer ==
    coveragePending => Active \/ RescanPending \/ pumpPending

EventuallyCompletesCoverage == <> (~coveragePending /\ state = "idle")

=======================================================================
