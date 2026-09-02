----------------------- MODULE ObolSubmissionPass -----------------------
\* One bounded model of the submission cursor, its full-rescan debt, and an
\* exact presentation-visibility delta.  Visibility is a planning input but
\* not an immutable mesh-inventory mutation.  A host synchronization which
\* observes a new visibility revision must nevertheless supply a level-
\* triggered producer until the exact delta is consumed.
\* Discovery may finish an exact append pass before compact inventory is
\* complete.  That pauses the cursor without consuming the rescan.  Fair
\* inventory publication and cursor ownership must eventually consume it.

EXTENDS TLC

States == {"idle", "active", "active-demand", "active-selective",
           "active-presentation", "active-allocation", "idle-rescan",
           "active-rescan"}

VARIABLES state, inventoryComplete, coveragePending, pumpPending,
          visibilityCurrent, presentationCurrent, allocationCurrent,
          capacityEvidenceCurrent, demandPending

vars == <<state, inventoryComplete, coveragePending, pumpPending,
          visibilityCurrent, presentationCurrent, allocationCurrent,
          capacityEvidenceCurrent, demandPending>>

Active == state \in {"active", "active-demand", "active-selective",
                     "active-presentation", "active-allocation",
                     "active-rescan"}
RescanPending == state \in {"idle-rescan", "active-rescan"}

Init ==
    /\ state \in {"idle", "idle-rescan", "active-rescan"}
    /\ inventoryComplete = (state = "idle")
    /\ coveragePending = TRUE
    /\ pumpPending = (state = "idle")
    /\ visibilityCurrent = TRUE
    /\ presentationCurrent = TRUE
    /\ allocationCurrent = TRUE
    /\ capacityEvidenceCurrent = TRUE
    /\ demandPending = FALSE

\* A retained hierarchy edit can change effective visibility without changing
\* asset identity or inventory.  synchronizePresentation is the durable host
\* boundary which observes the new revision and arms the progressive pump.
MutateVisibility ==
    /\ visibilityCurrent
    /\ visibilityCurrent' = FALSE
    /\ presentationCurrent' = FALSE
    /\ allocationCurrent' = FALSE
    /\ pumpPending' = TRUE
    /\ state' = IF state \in {"active-presentation", "active-allocation"}
                 THEN "active" ELSE state
    /\ UNCHANGED <<inventoryComplete, coveragePending,
                    capacityEvidenceCurrent, demandPending>>

PublishCompleteInventory ==
    /\ ~inventoryComplete
    /\ inventoryComplete' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<state, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

\* IDLE_RESCAN means the predecessor has already completed and the successor
\* full pass is owed.  Resuming therefore consumes the debt into ACTIVE;
\* ACTIVE_RESCAN instead means a still-running predecessor owes that successor.
ResumeRescan ==
    /\ state = "idle-rescan"
    /\ inventoryComplete
    /\ pumpPending
    /\ state' = "active"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

FinishRescanPass ==
    /\ state = "active-rescan"
    /\ state' = IF inventoryComplete THEN "active" ELSE "idle-rescan"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

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
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

FinishCoveragePass ==
    /\ state = "active"
    /\ coveragePending
    /\ state' = "idle"
    /\ coveragePending' = FALSE
    /\ pumpPending' =
        IF visibilityCurrent /\ presentationCurrent /\ allocationCurrent
        THEN FALSE ELSE TRUE
    /\ UNCHANGED <<inventoryComplete, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

\* A selective source/result pass may publish an immutable prefix which leaves
\* current-view quality unresolved.  It owns the shared cursor only until it
\* completes; the semantic demand obligation must survive that mechanical
\* owner and resume as one complete ordinary pass.
StartSelectivePass ==
    /\ state = "idle"
    /\ ~coveragePending
    /\ ~demandPending
    /\ visibilityCurrent
    /\ presentationCurrent
    /\ allocationCurrent
    /\ state' = "active-selective"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent, allocationCurrent,
                    capacityEvidenceCurrent, demandPending>>

PublishQualityDebt ==
    /\ state = "active-selective"
    /\ ~demandPending
    /\ demandPending' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<state, inventoryComplete, coveragePending,
                    visibilityCurrent, presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent>>

FinishSelectivePass ==
    /\ state = "active-selective"
    /\ state' = "idle"
    /\ pumpPending' = (demandPending \/ coveragePending \/
                        ~visibilityCurrent \/ ~presentationCurrent \/
                        ~allocationCurrent)
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent, allocationCurrent,
                    capacityEvidenceCurrent, demandPending>>

ResumeDemand ==
    /\ state = "idle"
    /\ demandPending
    /\ pumpPending
    /\ ~coveragePending
    /\ visibilityCurrent
    /\ presentationCurrent
    /\ allocationCurrent
    /\ state' = "active-demand"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent, allocationCurrent,
                    capacityEvidenceCurrent, demandPending>>

FinishDemandPass ==
    /\ state = "active-demand"
    /\ demandPending
    /\ state' = "idle"
    /\ demandPending' = FALSE
    /\ pumpPending' = (coveragePending \/ ~visibilityCurrent \/
                        ~presentationCurrent \/ ~allocationCurrent)
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent, allocationCurrent,
                    capacityEvidenceCurrent>>

ResumeVisibilityDelta ==
    /\ state = "idle"
    /\ ~visibilityCurrent
    /\ pumpPending
    /\ state' = "active"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

\* An exact delta may be consumed by a pass already active for coverage or by
\* the dedicated pass above.  It does not invalidate the complete denominator;
\* it updates only the changed entries and republishes that denominator.
FinishVisibilityDelta ==
    /\ Active
    /\ ~visibilityCurrent
    /\ visibilityCurrent' = TRUE
    /\ state' = IF state = "active" /\ ~coveragePending THEN "idle"
                 ELSE state
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<inventoryComplete, coveragePending, presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

\* The exact source census and the exact framebuffer are distinct facts.  The
\* source pass can restore an occurrence whose previous payload was retired;
\* the following exact frame must classify that occurrence as structural work
\* before allocation is allowed to choose its replacement cut.
ResumeVisibilityPresentation ==
    /\ state = "idle"
    /\ visibilityCurrent
    /\ ~presentationCurrent
    /\ pumpPending
    /\ state' = "active-presentation"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent, allocationCurrent,
                    capacityEvidenceCurrent, demandPending>>

FinishVisibilityPresentation ==
    /\ state = "active-presentation"
    /\ visibilityCurrent
    /\ ~presentationCurrent
    /\ state' = "idle"
    /\ presentationCurrent' = TRUE
    /\ pumpPending' = TRUE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

\* The selective source pass publishes effective visibility before minimax
\* examines the occurrence set.  Its successor reuses renderer-capacity
\* evidence but produces a new allocation certificate for that visibility
\* epoch.  Starting this action on MutateVisibility would inspect the old
\* hidden-instance set and omit newly restored entries.
ResumeVisibilityAllocation ==
    /\ state = "idle"
    /\ visibilityCurrent
    /\ presentationCurrent
    /\ ~allocationCurrent
    /\ pumpPending
    /\ state' = "active-allocation"
    /\ pumpPending' = FALSE
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    allocationCurrent, capacityEvidenceCurrent,
                    demandPending>>

FinishVisibilityAllocation ==
    /\ state = "active-allocation"
    /\ visibilityCurrent
    /\ presentationCurrent
    /\ ~allocationCurrent
    /\ state' = "idle"
    /\ allocationCurrent' = TRUE
    /\ pumpPending' = demandPending
    /\ UNCHANGED <<inventoryComplete, coveragePending, visibilityCurrent,
                    presentationCurrent,
                    capacityEvidenceCurrent, demandPending>>

Next == MutateVisibility \/ PublishCompleteInventory \/ ResumeRescan \/
        FinishRescanPass \/ ResumeCoverage \/ FinishCoveragePass \/
        StartSelectivePass \/ PublishQualityDebt \/ FinishSelectivePass \/
        ResumeDemand \/ FinishDemandPass \/
        ResumeVisibilityDelta \/ FinishVisibilityDelta \/
        ResumeVisibilityPresentation \/ FinishVisibilityPresentation \/
        ResumeVisibilityAllocation \/ FinishVisibilityAllocation

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishCompleteInventory)
    /\ WF_vars(ResumeRescan)
    /\ WF_vars(FinishRescanPass)
    /\ WF_vars(ResumeCoverage)
    /\ WF_vars(FinishCoveragePass)
    /\ WF_vars(PublishQualityDebt)
    /\ WF_vars(FinishSelectivePass)
    /\ WF_vars(ResumeDemand)
    /\ WF_vars(FinishDemandPass)
    /\ WF_vars(ResumeVisibilityDelta)
    /\ WF_vars(FinishVisibilityDelta)
    /\ WF_vars(ResumeVisibilityPresentation)
    /\ WF_vars(FinishVisibilityPresentation)
    /\ WF_vars(ResumeVisibilityAllocation)
    /\ WF_vars(FinishVisibilityAllocation)

TypeOK ==
    /\ state \in States
    /\ inventoryComplete \in BOOLEAN
    /\ coveragePending \in BOOLEAN
    /\ pumpPending \in BOOLEAN
    /\ visibilityCurrent \in BOOLEAN
    /\ presentationCurrent \in BOOLEAN
    /\ allocationCurrent \in BOOLEAN
    /\ capacityEvidenceCurrent \in BOOLEAN
    /\ demandPending \in BOOLEAN

PausedPassPreservesDebt == state = "idle-rescan" => ~Active /\ RescanPending

ActiveRescanOwnsBothFacts ==
    state = "active-rescan" => Active /\ RescanPending

EventuallyConsumesRescan == <> (state = "idle" /\ ~RescanPending)

CoverageHasProducer ==
    coveragePending => Active \/ RescanPending \/ pumpPending

DemandHasProducer ==
    demandPending => Active \/ pumpPending

VisibilityDeltaHasProducer ==
    ~visibilityCurrent => Active \/ RescanPending \/ pumpPending

VisibilityAllocationHasProducer ==
    visibilityCurrent /\ presentationCurrent /\ ~allocationCurrent =>
        Active \/ RescanPending \/ pumpPending

VisibilityPresentationHasProducer ==
    visibilityCurrent /\ ~presentationCurrent =>
        Active \/ RescanPending \/ pumpPending

AllocationRequiresExactVisibilityPresentation ==
    allocationCurrent => visibilityCurrent /\ presentationCurrent

VisibilityPreservesCapacityEvidence == capacityEvidenceCurrent

Terminal == state = "idle" /\ ~Active
DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyCompletesCoverage == <> (~coveragePending /\ state = "idle")

\* This local model admits unbounded external visibility mutations.  A pending
\* demand pass must therefore either complete or be superseded by a stronger
\* visibility/presentation/allocation epoch; the composed convergence model
\* separately proves termination after external input closes.
DemandEventuallyConsumesOrIsSuperseded ==
    [](demandPending ~>
        (~demandPending \/ ~visibilityCurrent \/ ~presentationCurrent \/
         ~allocationCurrent))

EventuallyAppliesVisibilityDelta ==
    [](~visibilityCurrent ~> visibilityCurrent)

QuiescentVisibilityEventuallyPresented ==
    (<>[]visibilityCurrent) => (<>[]presentationCurrent)

\* No renderer can promise a terminal allocation while the application keeps
\* replacing its visibility input forever.  Once that input becomes quiet,
\* the last exact delta and its successor allocation must both become stable.
QuiescentVisibilityEventuallyReallocated ==
    (<>[]visibilityCurrent) => (<>[](presentationCurrent /\ allocationCurrent))

=======================================================================
