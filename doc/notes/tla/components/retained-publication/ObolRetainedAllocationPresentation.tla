---------------- MODULE ObolRetainedAllocationPresentation ----------------
\* A retained allocation can mutate three independent presentation channels:
\* occurrence-local PoP cuts, spatial page/channel bindings, and the
\* point/mesh protection set consumed by Obol's camera-local classifier.
\* Committing any of them is not evidence that it reached a framebuffer.  In
\* particular, a numerically unchanged cut may still lack newly visible pages,
\* and point-only allocations have no changed CadPayload cut from which the
\* ordinary render barrier could infer their obligation.
\*
\* This model isolates that handoff.  A changed protection set must arm one
\* frame producer, complete classification for the committed generation, and
\* execute an exact frame before the allocation may be called realized.  The
\* retained work outside occurrence-local CAD allocation is an environmental
\* input: applying either presentation channel must not derive or rewrite it
\* from the resulting framebuffer.

EXTENDS Naturals, TLC

Phases == {"idle", "committed", "presented", "realized"}
AbstractCostMaximum == 2

VARIABLES phase, protectionChanged, presentationRequested, cutsApplied,
          spatialPresentationReady, classifierCurrent, exactFrame,
          unmanagedPresentationCost, allocationExternalCost

vars == <<phase, protectionChanged, presentationRequested, cutsApplied,
          spatialPresentationReady, classifierCurrent, exactFrame,
          unmanagedPresentationCost, allocationExternalCost>>

costs == <<unmanagedPresentationCost, allocationExternalCost>>

Init ==
    /\ phase = "idle"
    /\ protectionChanged \in BOOLEAN
    /\ presentationRequested = FALSE
    /\ cutsApplied = FALSE
    /\ spatialPresentationReady \in BOOLEAN
    /\ classifierCurrent = TRUE
    /\ exactFrame \in BOOLEAN
    /\ unmanagedPresentationCost \in 0..AbstractCostMaximum
    /\ allocationExternalCost = unmanagedPresentationCost

CommitAllocation ==
    /\ phase = "idle"
    /\ phase' = "committed"
    /\ presentationRequested' =
        (protectionChanged \/ ~spatialPresentationReady \/ ~exactFrame)
    /\ cutsApplied' = TRUE
    /\ classifierCurrent' = ~protectionChanged
    /\ exactFrame' =
        (exactFrame /\ ~protectionChanged /\ spatialPresentationReady)
    /\ UNCHANGED <<protectionChanged, spatialPresentationReady, costs>>

CompleteSpatialPresentation ==
    /\ phase = "committed"
    /\ presentationRequested
    /\ ~spatialPresentationReady
    /\ spatialPresentationReady' = TRUE
    /\ UNCHANGED <<phase, protectionChanged, presentationRequested,
                    cutsApplied, classifierCurrent, exactFrame,
                    unmanagedPresentationCost, allocationExternalCost>>

CompleteClassification ==
    /\ phase = "committed"
    /\ protectionChanged
    /\ presentationRequested
    /\ classifierCurrent' = TRUE
    /\ UNCHANGED <<phase, protectionChanged, presentationRequested,
                    cutsApplied, spatialPresentationReady, exactFrame,
                    unmanagedPresentationCost,
                    allocationExternalCost>>

PresentExactFrame ==
    /\ phase = "committed"
    /\ presentationRequested
    /\ classifierCurrent
    /\ spatialPresentationReady
    /\ phase' = "presented"
    /\ exactFrame' = TRUE
    /\ UNCHANGED <<protectionChanged, presentationRequested, cutsApplied,
                    spatialPresentationReady, classifierCurrent,
                    unmanagedPresentationCost,
                    allocationExternalCost>>

RealizeUnchangedAllocation ==
    /\ phase = "committed"
    /\ ~presentationRequested
    /\ cutsApplied
    /\ spatialPresentationReady
    /\ classifierCurrent
    /\ exactFrame
    /\ phase' = "realized"
    /\ UNCHANGED <<protectionChanged, presentationRequested, cutsApplied,
                    spatialPresentationReady, classifierCurrent, exactFrame,
                    unmanagedPresentationCost, allocationExternalCost>>

RealizePresentedAllocation ==
    /\ phase = "presented"
    /\ cutsApplied
    /\ spatialPresentationReady
    /\ classifierCurrent
    /\ exactFrame
    /\ phase' = "realized"
    /\ presentationRequested' = FALSE
    /\ UNCHANGED <<protectionChanged, cutsApplied,
                    spatialPresentationReady, classifierCurrent, exactFrame,
                    unmanagedPresentationCost,
                    allocationExternalCost>>

Next == CommitAllocation \/ CompleteSpatialPresentation \/
        CompleteClassification \/ PresentExactFrame \/
        RealizeUnchangedAllocation \/ RealizePresentedAllocation

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(CommitAllocation)
    /\ WF_vars(CompleteSpatialPresentation)
    /\ WF_vars(CompleteClassification)
    /\ WF_vars(PresentExactFrame)
    /\ WF_vars(RealizeUnchangedAllocation)
    /\ WF_vars(RealizePresentedAllocation)

TypeOK ==
    /\ phase \in Phases
    /\ protectionChanged \in BOOLEAN
    /\ presentationRequested \in BOOLEAN
    /\ cutsApplied \in BOOLEAN
    /\ spatialPresentationReady \in BOOLEAN
    /\ classifierCurrent \in BOOLEAN
    /\ exactFrame \in BOOLEAN
    /\ unmanagedPresentationCost \in 0..AbstractCostMaximum
    /\ allocationExternalCost \in 0..AbstractCostMaximum

ChangedCommitHasProducer ==
    phase = "committed" /\
    (~spatialPresentationReady \/ ~classifierCurrent \/ ~exactFrame) =>
        presentationRequested

RealizedAllocationWasPresented ==
    phase = "realized" =>
        cutsApplied /\ spatialPresentationReady /\ classifierCurrent /\
        exactFrame

ExternalAllocationInputStable ==
    allocationExternalCost = unmanagedPresentationCost

Terminal == phase = "realized"
DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyRealized == <> (phase = "realized")

=============================================================================
