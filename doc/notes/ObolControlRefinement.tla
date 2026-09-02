---------------------- MODULE ObolControlRefinement ----------------------
\* Executable control-refinement contract shared by the production C++ map
\* and ObolProgressivePipeline.tla.  Every production input is represented
\* independently here and refines to one of nine finite obligations.  Keeping
\* the map one-to-one is intentional: adding an implementation fact without
\* updating this model must be visible as contract drift even when two facts
\* currently refine to the same obligation.  This model deliberately excludes
\* numeric planning and renderer work.

EXTENDS Naturals, FiniteSets, TLC

Facts == {
    "interaction", "inventory", "availability", "result", "publication",
    "submission", "demandRefresh", "submissionRescan", "submissionDelta",
    "qualityProbe",
    "retainedAllocation", "retainedAllocationTransaction",
    "importanceCensus", "residentAdmissionRetry", "capacityAllocation",
    "residentGrowth", "pointTriangleRecovery", "structuralFrontier",
    "presentationReplay", "exactPresentation", "presentationBarrier",
    "capacityFrame",
    "pointAdmissionFrame", "pointCalibration", "capacityCalibration",
    "headroomProbe", "handoff", "compaction", "cacheWrite"
}

Work == {
    "interaction", "inventory", "availability", "publication", "planning",
    "presentation", "handoff", "compaction", "cacheWrite"
}

Background == {"compaction", "cacheWrite"}

PlanningFacts == {
    "submission", "demandRefresh", "submissionRescan", "submissionDelta",
    "qualityProbe",
    "retainedAllocation", "retainedAllocationTransaction",
    "importanceCensus", "residentAdmissionRetry", "capacityAllocation",
    "residentGrowth", "pointTriangleRecovery", "structuralFrontier"
}

PresentationFacts == {
    "presentationReplay", "exactPresentation", "presentationBarrier",
    "capacityFrame",
    "pointAdmissionFrame", "pointCalibration", "capacityCalibration",
    "headroomProbe"
}

\* Aliases inside PlanningFacts and PresentationFacts cannot change owner
\* precedence.  Explore all combinations of the semantically distinct
\* classes plus presentationReplay, whose closed transaction deliberately
\* precedes every other owner, and capacityAllocation, whose exact frame is a
\* downstream successor.  CorrectFactMap below still checks every concrete
\* production fact independently.
RepresentativeFacts == {
    "interaction", "inventory", "availability", "result", "submission",
    "capacityAllocation", "presentationReplay", "exactPresentation",
    "presentationBarrier", "handoff", "compaction", "cacheWrite"
}

FactWork(fact) ==
    IF fact = "interaction" THEN "interaction"
    ELSE IF fact = "inventory" THEN "inventory"
    ELSE IF fact = "availability" THEN "availability"
    ELSE IF fact \in {"result", "publication"} THEN "publication"
    ELSE IF fact \in PlanningFacts THEN "planning"
    ELSE IF fact \in PresentationFacts THEN "presentation"
    ELSE IF fact = "handoff" THEN "handoff"
    ELSE IF fact = "compaction" THEN "compaction"
    ELSE "cacheWrite"

Obligations(factSet) == {FactWork(fact) : fact \in factSet}

Owner(factSet) ==
    IF "presentationReplay" \in factSet THEN "presentation"
    ELSE IF "interaction" \in factSet THEN "interaction"
    ELSE IF factSet \intersect {"result", "publication"} # {} THEN
        "publication"
    ELSE IF "inventory" \in factSet THEN "inventory"
    ELSE IF "availability" \in factSet THEN "availability"
    ELSE IF "capacityAllocation" \in factSet THEN "planning"
    ELSE IF "presentation" \in Obligations(factSet) THEN "presentation"
    ELSE IF "planning" \in Obligations(factSet) THEN "planning"
    ELSE IF "handoff" \in factSet THEN "handoff"
    ELSE IF "compaction" \in factSet THEN "compaction"
    ELSE IF "cacheWrite" \in factSet THEN "cacheWrite"
    ELSE "none"

VARIABLES facts, error, probeFact

vars == <<facts, error, probeFact>>

Foreground == Obligations(facts) \ Background
Terminal == Foreground = {}
Ready == Terminal /\ ~error

TypeOK ==
    /\ facts \subseteq RepresentativeFacts
    /\ error \in BOOLEAN
    /\ probeFact \in Facts
    /\ Obligations(facts) \subseteq Work

CorrectFactMap ==
    /\ FactWork("interaction") = "interaction"
    /\ FactWork("inventory") = "inventory"
    /\ FactWork("availability") = "availability"
    /\ \A fact \in {"result", "publication"}:
           FactWork(fact) = "publication"
    /\ \A fact \in PlanningFacts: FactWork(fact) = "planning"
    /\ \A fact \in PresentationFacts: FactWork(fact) = "presentation"
    /\ FactWork("handoff") = "handoff"
    /\ FactWork("compaction") = "compaction"
    /\ FactWork("cacheWrite") = "cacheWrite"
    /\ FactWork(probeFact) \in Work

OwnerIsValid ==
    /\ (Owner(facts) = "none") = (facts = {})
    /\ Owner(facts) # "none" => Owner(facts) \in Obligations(facts)

TerminalHasNoForegroundWork == Terminal => Foreground = {}
ReadyIsValid == Ready => Terminal /\ ~error

Init ==
    /\ facts \in SUBSET RepresentativeFacts
    /\ error \in BOOLEAN
    /\ probeFact \in Facts

DischargeOwner ==
    /\ Owner(facts) # "none"
    /\ facts' = {fact \in facts : FactWork(fact) # Owner(facts)}
    /\ UNCHANGED <<error, probeFact>>

Next == DischargeOwner

Spec == Init /\ [][Next]_vars /\ WF_vars(DischargeOwner)

DischargeStrictlyReduces ==
    [][DischargeOwner => Cardinality(facts') < Cardinality(facts)]_vars

EventuallyAllWorkRetires == <> (facts = {})

=============================================================================
