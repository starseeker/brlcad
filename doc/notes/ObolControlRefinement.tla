---------------------- MODULE ObolControlRefinement ----------------------
\* Executable control-refinement contract shared by the production C++ map
\* and ObolProgressivePipeline.tla.  Twenty fact classes represent the
\* twenty-one remaining implementation inputs and refine to nine finite
\* obligations.  The production structuralFrontier and pointTriangleRecovery
\* booleans are symmetry-equivalent here: both map only to planning, while the
\* exhaustive C++ refinement test retains and checks them independently.  This
\* model deliberately excludes numeric planning and renderer work.

EXTENDS Naturals, FiniteSets, TLC

Facts == {
    "interaction", "inventory", "availability", "result", "publication",
    "submission", "submissionRescan", "retainedAllocation",
    "residentGrowth", "pointTriangleRecovery", "presentationReplay",
    "presentationBarrier", "capacityFrame", "pointAdmissionFrame",
    "pointCalibration", "capacityCalibration", "headroomProbe", "handoff",
    "compaction", "cacheWrite"
}

Work == {
    "interaction", "inventory", "availability", "publication", "planning",
    "presentation", "handoff", "compaction", "cacheWrite"
}

Background == {"compaction", "cacheWrite"}

FactWork(fact) ==
    IF fact = "interaction" THEN "interaction"
    ELSE IF fact = "inventory" THEN "inventory"
    ELSE IF fact = "availability" THEN "availability"
    ELSE IF fact \in {"result", "publication"} THEN "publication"
    ELSE IF fact \in {"submission", "submissionRescan",
                       "retainedAllocation", "residentGrowth",
                       "pointTriangleRecovery"} THEN "planning"
    ELSE IF fact \in {"presentationReplay", "presentationBarrier",
                       "capacityFrame", "pointAdmissionFrame",
                       "pointCalibration", "capacityCalibration",
                       "headroomProbe"} THEN "presentation"
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
    ELSE IF "presentation" \in Obligations(factSet) THEN "presentation"
    ELSE IF "planning" \in Obligations(factSet) THEN "planning"
    ELSE IF "handoff" \in factSet THEN "handoff"
    ELSE IF "compaction" \in factSet THEN "compaction"
    ELSE IF "cacheWrite" \in factSet THEN "cacheWrite"
    ELSE "none"

VARIABLES facts, error

vars == <<facts, error>>

Foreground == Obligations(facts) \ Background
Terminal == Foreground = {}
Ready == Terminal /\ ~error

TypeOK ==
    /\ facts \subseteq Facts
    /\ error \in BOOLEAN
    /\ Obligations(facts) \subseteq Work

OwnerIsValid ==
    /\ (Owner(facts) = "none") = (facts = {})
    /\ Owner(facts) # "none" => Owner(facts) \in Obligations(facts)

TerminalHasNoForegroundWork == Terminal => Foreground = {}
ReadyIsValid == Ready => Terminal /\ ~error

Init ==
    /\ facts \in SUBSET Facts
    /\ error \in BOOLEAN

DischargeOwner ==
    /\ Owner(facts) # "none"
    /\ facts' = {fact \in facts : FactWork(fact) # Owner(facts)}
    /\ UNCHANGED error

Next == DischargeOwner

Spec == Init /\ [][Next]_vars /\ WF_vars(DischargeOwner)

DischargeStrictlyReduces ==
    [][DischargeOwner => Cardinality(facts') < Cardinality(facts)]_vars

EventuallyAllWorkRetires == <> (facts = {})

=============================================================================
