---------------- MODULE ObolAssetPublicationComposition ----------------
\* Composition seam between one demand-independent immutable asset producer,
\* demand-qualified live publication, and durable cache completion.
\*
\* A camera/policy demand may overtake any live page without cancelling the
\* expensive asset build.  Ready pages remain immutable producer output; the
\* owner thread binds them to the newest demand immediately before publication.
\* A stale page is useful only as an intermediate retained layer and cannot
\* satisfy current-demand completion.  Whole-object coverage remains visible
\* until the complete hierarchy is available.  Durable cache discovery is
\* permitted only after every page is complete, independently of which pages
\* happened to be presented during cold construction.

EXTENDS Naturals, FiniteSets, TLC

CONSTANT MaxDemand

Pages == {"first", "second"}
Demands == 0..MaxDemand
NoneDemand == MaxDemand + 1
OptionalDemand == 0..NoneDemand
ProducerStates ==
    {"building", "complete", "constrained", "failed", "cancelled"}

VARIABLES inputOpen,
          demand,
          producer,
          readyPages,
          pageDemand,
          finalGeometry,
          finalPresentedDemand,
          coverage,
          cacheMarked,
          constraintDemand,
          failureDemand,
          framePending

vars == <<inputOpen, demand, producer, readyPages, pageDemand,
          finalGeometry, finalPresentedDemand, coverage, cacheMarked,
          constraintDemand, failureDemand, framePending>>

TypeOK ==
    /\ inputOpen \in BOOLEAN
    /\ demand \in Demands
    /\ producer \in ProducerStates
    /\ readyPages \subseteq Pages
    /\ pageDemand \in [Pages -> OptionalDemand]
    /\ finalGeometry \in BOOLEAN
    /\ finalPresentedDemand \in OptionalDemand
    /\ coverage \in BOOLEAN
    /\ cacheMarked \in BOOLEAN
    /\ constraintDemand \in OptionalDemand
    /\ failureDemand \in OptionalDemand
    /\ framePending \in BOOLEAN

Init ==
    /\ inputOpen = TRUE
    /\ demand = 0
    /\ producer = "building"
    /\ readyPages = {}
    /\ pageDemand = [page \in Pages |-> NoneDemand]
    /\ finalGeometry = FALSE
    /\ finalPresentedDemand = NoneDemand
    /\ coverage = TRUE
    /\ cacheMarked = FALSE
    /\ constraintDemand = NoneDemand
    /\ failureDemand = NoneDemand
    /\ framePending = FALSE

MoveDemand ==
    /\ inputOpen
    /\ demand < MaxDemand
    /\ demand' = demand + 1
    /\ constraintDemand' = NoneDemand
    /\ failureDemand' = NoneDemand
    /\ UNCHANGED <<inputOpen, producer, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, framePending>>

\* A terminal outcome belongs to the request which observed it, not to
\* immutable pages already completed by the asset producer.  Once a newer
\* demand has retired that outcome, the same producer resumes from its
\* retained prefix.
RestartConstrained ==
    /\ producer = "constrained"
    /\ constraintDemand # demand
    /\ producer' = "building"
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand,
                    framePending>>

RestartFailed ==
    /\ producer = "failed"
    /\ failureDemand # demand
    /\ producer' = "building"
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand,
                    framePending>>

CloseInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ UNCHANGED <<demand, producer, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand,
                    framePending>>

CompletePage(page) ==
    /\ producer = "building"
    /\ page \in Pages \ readyPages
    /\ readyPages' = readyPages \cup {page}
    /\ UNCHANGED <<inputOpen, demand, producer, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand,
                    framePending>>

\* Binding is a current owner-thread act.  It does not mutate immutable page
\* bytes and an old binding cannot become a failure for a newer demand.
PublishPage(page) ==
    /\ page \in readyPages
    /\ pageDemand[page] # demand
    /\ pageDemand' = [pageDemand EXCEPT ![page] = demand]
    /\ framePending' = TRUE
    /\ UNCHANGED <<inputOpen, demand, producer, readyPages,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand>>

FinalizeAsset ==
    /\ producer = "building"
    /\ readyPages = Pages
    /\ producer' = "complete"
    /\ finalGeometry' = TRUE
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalPresentedDemand, coverage, cacheMarked,
                    constraintDemand, failureDemand, framePending>>

\* The complete hierarchy is globally representative; binding it to the
\* newest demand may retire temporary whole-object coverage without replaying
\* every live-page publication from earlier demands.
PublishFinal ==
    /\ finalGeometry
    /\ finalPresentedDemand # demand
    /\ finalPresentedDemand' = demand
    /\ coverage' = FALSE
    /\ framePending' = TRUE
    /\ UNCHANGED <<inputOpen, demand, producer, readyPages, pageDemand,
                    finalGeometry, cacheMarked, constraintDemand,
                    failureDemand>>

MarkCache ==
    /\ producer = "complete"
    /\ finalGeometry
    /\ readyPages = Pages
    /\ ~cacheMarked
    /\ cacheMarked' = TRUE
    /\ UNCHANGED <<inputOpen, demand, producer, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    constraintDemand, failureDemand, framePending>>

Constrain ==
    /\ producer = "building"
    /\ producer' = "constrained"
    /\ constraintDemand' = demand
    /\ framePending' = TRUE
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, failureDemand>>

\* Invalid geometry or an unrecoverable provider error is not a quality or
\* resource constraint.  It carries independent current-demand evidence and
\* reaches Failed rather than being relabeled Constrained.
Fail ==
    /\ producer = "building"
    /\ producer' = "failed"
    /\ failureDemand' = demand
    /\ framePending' = TRUE
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand>>

Cancel ==
    /\ producer = "building"
    /\ producer' = "cancelled"
    /\ UNCHANGED <<inputOpen, demand, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand,
                    framePending>>

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ UNCHANGED <<inputOpen, demand, producer, readyPages, pageDemand,
                    finalGeometry, finalPresentedDemand, coverage,
                    cacheMarked, constraintDemand, failureDemand>>

Quiescent ==
    /\ ~inputOpen
    /\ ~framePending
    /\ producer \in {"complete", "constrained", "failed", "cancelled"}
    /\ UNCHANGED vars

Next ==
    \/ MoveDemand
    \/ RestartConstrained
    \/ RestartFailed
    \/ CloseInput
    \/ \E page \in Pages: CompletePage(page)
    \/ \E page \in Pages: PublishPage(page)
    \/ FinalizeAsset
    \/ PublishFinal
    \/ MarkCache
    \/ Constrain
    \/ Fail
    \/ Cancel
    \/ CompleteFrame
    \/ Quiescent

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(RestartConstrained)
    /\ WF_vars(RestartFailed)
    /\ \A page \in Pages: WF_vars(CompletePage(page))
    /\ WF_vars(FinalizeAsset)
    /\ WF_vars(PublishFinal)
    /\ WF_vars(MarkCache)
    /\ WF_vars(CompleteFrame)

ReadyBeforePublication ==
    \A page \in Pages:
        pageDemand[page] # NoneDemand => page \in readyPages

CacheRequiresCompleteAsset ==
    cacheMarked =>
        /\ producer = "complete"
        /\ finalGeometry
        /\ readyPages = Pages

CoverageRetiresOnlyForFinalGeometry == ~coverage => finalGeometry

BuildingHasNoTerminalArtifact ==
    producer = "building" => ~finalGeometry /\ ~cacheMarked

CurrentConstraintOnly ==
    constraintDemand = NoneDemand \/ constraintDemand = demand

CurrentFailureOnly == failureDemand = NoneDemand \/ failureDemand = demand

TypedOutcomeIsExclusive ==
    ~(constraintDemand = demand /\ failureDemand = demand)

AlwaysDisplayable == coverage \/ finalGeometry \/
    (\E page \in Pages: pageDemand[page] # NoneDemand)

Terminal ==
    /\ ~inputOpen
    /\ ~framePending
    /\ producer \in {"complete", "constrained", "failed", "cancelled"}
    /\ (finalPresentedDemand = demand \/ constraintDemand = demand \/
         failureDemand = demand \/
         producer = "cancelled")

\* Cancellation is a successful terminal control outcome only after the
\* input scope has closed and therefore owns no current publication demand.
Ready ==
    /\ Terminal
    /\ (finalPresentedDemand = demand \/ producer = "cancelled")
Constrained == Terminal /\ constraintDemand = demand
Failed == Terminal /\ failureDemand = demand

EventuallyTerminalAfterInputCloses ==
    [](~inputOpen => <>Terminal)

=============================================================================
