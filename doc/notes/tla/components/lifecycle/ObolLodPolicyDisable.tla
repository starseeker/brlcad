----------------------- MODULE ObolLodPolicyDisable -----------------------
\* Focused ownership contract for toggling automatic mesh LoD.
\*
\* A policy-off transition retires every automatic owner atomically.  It does
\* not delete the immutable presentation already on screen, cancel an
\* independent database provider, or suppress the presentation-only repaint
\* needed to show the policy change.  Camera bookkeeping remains active while
\* policy is off, but it may not manufacture automatic demand.
\*
\* AutomaticDomains is the stable reducer-level map, not a second phase
\* machine.  The production fields represented by each domain are documented
\* in libbobol_formal_models.md and checked independently by the C++ retirement
\* test.  Geometry and numeric quality are deliberately abstracted.

EXTENDS Naturals, FiniteSets, TLC

CONSTANT MaxPolicyEpoch

AutomaticDomains == {
    "service", "availability", "submission", "planning", "viewDemand",
    "coverage", "interaction", "presentation", "pointQuality",
    "staticQuality", "structuralRepair", "compaction", "capacityHost"
}

VARIABLES automatic,
          policyEpoch,
          retainedPresentation,
          automaticWork,
          repaintPending,
          providerPending

vars == <<automatic, policyEpoch, retainedPresentation, automaticWork,
          repaintPending, providerPending>>

TypeOK ==
    /\ automatic \in BOOLEAN
    /\ policyEpoch \in 0..MaxPolicyEpoch
    /\ retainedPresentation \in BOOLEAN
    /\ automaticWork \subseteq AutomaticDomains
    /\ repaintPending \in BOOLEAN
    /\ providerPending \in BOOLEAN

Init ==
    /\ automatic = TRUE
    /\ policyEpoch = 0
    /\ retainedPresentation = TRUE
    /\ automaticWork \in SUBSET AutomaticDomains
    /\ repaintPending \in BOOLEAN
    /\ providerPending \in BOOLEAN

\* Odd transitions disable; even transitions re-enable.  Disable is one
\* reducer transition.  Re-enable starts only the minimum fresh-view work;
\* every other domain must be armed by its owning producer.
TogglePolicy ==
    /\ policyEpoch < MaxPolicyEpoch
    /\ policyEpoch' = policyEpoch + 1
    /\ automatic' = ~automatic
    /\ retainedPresentation' = retainedPresentation
    /\ repaintPending' = TRUE
    /\ providerPending' = providerPending
    /\ automaticWork' =
        IF automatic
        THEN {}
        ELSE {"service", "submission", "coverage", "capacityHost"}

\* Automatic policy may accumulate any combination of reducer domains before
\* a later toggle.  Modeling them as a set makes the atomic retirement rule
\* explicit and keeps adding an implementation owner from silently escaping
\* the policy-off transaction.
ArmAutomaticWork(domain) ==
    /\ automatic
    /\ domain \in AutomaticDomains
    /\ domain \notin automaticWork
    /\ automaticWork' = automaticWork \union {domain}
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    repaintPending, providerPending>>

CompleteAutomaticWork(domain) ==
    /\ automatic
    /\ domain \in automaticWork
    /\ automaticWork' = automaticWork \ {domain}
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    repaintPending, providerPending>>

\* A view revision is still recorded while policy is disabled.  Only enabled
\* policy may translate it into view-demand work.  This is the formal guard
\* for advanceLodViewRevision and similar camera-side entry points.
CameraChange ==
    /\ ~repaintPending \/
       (automatic /\ "viewDemand" \notin automaticWork)
    /\ repaintPending' = TRUE
    /\ automaticWork' =
        IF automatic
        THEN automaticWork \union {"viewDemand"}
        ELSE automaticWork
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    providerPending>>

CompleteProvider ==
    /\ providerPending
    /\ providerPending' = FALSE
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    automaticWork, repaintPending>>

PresentRepaint ==
    /\ repaintPending
    /\ repaintPending' = FALSE
    /\ UNCHANGED <<automatic, policyEpoch, retainedPresentation,
                    automaticWork, providerPending>>

Next ==
    \/ TogglePolicy
    \/ \E domain \in AutomaticDomains: ArmAutomaticWork(domain)
    \/ \E domain \in AutomaticDomains: CompleteAutomaticWork(domain)
    \/ CameraChange
    \/ CompleteProvider
    \/ PresentRepaint

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ \A domain \in AutomaticDomains:
           WF_vars(CompleteAutomaticWork(domain))
    /\ WF_vars(CompleteProvider)
    /\ WF_vars(PresentRepaint)

NoAutomaticOwner == automaticWork = {}

RetainedPresentationSurvives == retainedPresentation

DisabledOwnsNoAutomaticWork == ~automatic => NoAutomaticOwner

DisabledRenderIsPresentationOnly ==
    ~automatic => "capacityHost" \notin automaticWork

DisabledCameraCreatesNoDemand ==
    ~automatic => "viewDemand" \notin automaticWork

Terminal ==
    /\ policyEpoch = MaxPolicyEpoch
    /\ ~automatic
    /\ NoAutomaticOwner
    /\ ~providerPending
    /\ ~repaintPending

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyQuiescentAfterFinalDisable ==
    []((policyEpoch = MaxPolicyEpoch /\ ~automatic) =>
       <> (NoAutomaticOwner /\ ~providerPending /\ ~repaintPending))

=============================================================================
