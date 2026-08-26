----------------------- MODULE ObolSubmissionPass -----------------------
\* One bounded model of the submission cursor and its full-rescan debt.
\* Discovery may finish an exact append pass before compact inventory is
\* complete.  That pauses the cursor without consuming the rescan.  Fair
\* inventory publication and cursor ownership must eventually consume it.

EXTENDS TLC

States == {"idle", "active", "idle-rescan", "active-rescan"}

VARIABLES state, inventoryComplete

vars == <<state, inventoryComplete>>

Active == state \in {"active", "active-rescan"}
RescanPending == state \in {"idle-rescan", "active-rescan"}

Init ==
    /\ state = "idle-rescan"
    /\ inventoryComplete = FALSE

PublishCompleteInventory ==
    /\ ~inventoryComplete
    /\ inventoryComplete' = TRUE
    /\ UNCHANGED state

ResumeRescan ==
    /\ state = "idle-rescan"
    /\ state' = "active-rescan"
    /\ UNCHANGED inventoryComplete

FinishRescanPass ==
    /\ state = "active-rescan"
    /\ state' = IF inventoryComplete THEN "idle" ELSE "idle-rescan"
    /\ UNCHANGED inventoryComplete

Next == PublishCompleteInventory \/ ResumeRescan \/ FinishRescanPass

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishCompleteInventory)
    /\ WF_vars(ResumeRescan)
    /\ WF_vars(FinishRescanPass)

TypeOK == state \in States /\ inventoryComplete \in BOOLEAN

PausedPassPreservesDebt == state = "idle-rescan" => ~Active /\ RescanPending

ActiveRescanOwnsBothFacts ==
    state = "active-rescan" => Active /\ RescanPending

EventuallyConsumesRescan == <> (state = "idle" /\ ~RescanPending)

=======================================================================
