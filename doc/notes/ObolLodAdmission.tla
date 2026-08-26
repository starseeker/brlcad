--------------------------- MODULE ObolLodAdmission -------------------------
\* Bounded admission model for the compact source -> view allocator boundary.
\*
\* A source profile is deliberately a conservative fast-path gate, not a
\* quality decision.  Safe scenes may admit materially visible small meshes
\* directly; large scenes use the normal progressive path.  In both cases a
\* subpixel occurrence is represented by the aggregate channel and consumes
\* neither mesh draw budget nor source-preparation capacity.
\*
\* This model does not represent rendering, I/O, or PoP cuts.  It verifies the
\* policy combinations that must remain true before those implementations run:
\* no direct escape from a large scene, no mesh demand for a subpixel leaf,
\* bounded visible admission, and a terminal constrained state rather than an
\* allocation/reclassification cycle.

EXTENDS Naturals, TLC

Profiles == {"safe", "large"}
Candidates == {"hero", "ordinary", "tiny"}
Tiny == "tiny"
Significant == Candidates \ {Tiny}
Modes == {"box", "aggregate", "direct", "pop"}
InitialBudgets == 0..3

Cost(candidate) == IF candidate = "hero" THEN 2 ELSE 1

VARIABLES profile,
          budget,
          mode,
          allocationPending,
          constrained

vars == <<profile, budget, mode, allocationPending, constrained>>

TypeOK ==
    /\ profile \in Profiles
    /\ budget \in 0..3
    /\ mode \in [Candidates -> Modes]
    /\ allocationPending \in BOOLEAN
    /\ constrained \in BOOLEAN

Init ==
    /\ profile \in Profiles
    /\ budget \in InitialBudgets
    /\ mode = [candidate \in Candidates |-> "box"]
    /\ allocationPending = TRUE
    /\ constrained = FALSE

ClassifyTiny ==
    /\ allocationPending
    /\ mode[Tiny] = "box"
    /\ mode' = [mode EXCEPT ![Tiny] = "aggregate"]
    /\ UNCHANGED <<profile, budget, allocationPending, constrained>>

AdmitDirect(candidate) ==
    /\ candidate \in Significant
    /\ profile = "safe"
    /\ allocationPending
    /\ mode[candidate] = "box"
    /\ budget >= Cost(candidate)
    /\ candidate = "hero" \/ mode["hero"] # "box" \/
       budget < Cost("hero")
    /\ mode' = [mode EXCEPT ![candidate] = "direct"]
    /\ budget' = budget - Cost(candidate)
    /\ UNCHANGED <<profile, allocationPending, constrained>>

AdmitPop(candidate) ==
    /\ candidate \in Significant
    /\ profile = "large"
    /\ allocationPending
    /\ mode[candidate] = "box"
    /\ budget >= Cost(candidate)
    /\ candidate = "hero" \/ mode["hero"] # "box" \/
       budget < Cost("hero")
    /\ mode' = [mode EXCEPT ![candidate] = "pop"]
    /\ budget' = budget - Cost(candidate)
    /\ UNCHANGED <<profile, allocationPending, constrained>>

AllClassified == \A candidate \in Candidates: mode[candidate] # "box"
AnyAffordable ==
    \E candidate \in Significant:
        /\ mode[candidate] = "box"
        /\ budget >= Cost(candidate)

Complete ==
    /\ allocationPending
    /\ AllClassified
    /\ allocationPending' = FALSE
    /\ UNCHANGED <<profile, budget, mode, constrained>>

Constrain ==
    /\ allocationPending
    /\ mode[Tiny] # "box"
    /\ ~AnyAffordable
    /\ ~AllClassified
    /\ allocationPending' = FALSE
    /\ constrained' = TRUE
    /\ UNCHANGED <<profile, budget, mode>>

Next ==
    \/ ClassifyTiny
    \/ \E candidate \in Significant: AdmitDirect(candidate)
    \/ \E candidate \in Significant: AdmitPop(candidate)
    \/ Complete
    \/ Constrain

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ClassifyTiny)
    /\ \A candidate \in Significant: WF_vars(AdmitDirect(candidate))
    /\ \A candidate \in Significant: WF_vars(AdmitPop(candidate))
    /\ WF_vars(Complete)
    /\ WF_vars(Constrain)

BudgetSafety == budget >= 0

\* A large scene never obtains a source-size-only direct escape hatch.
LargeSceneHasNoDirectMeshes ==
    profile = "large" =>
        \A candidate \in Candidates: mode[candidate] # "direct"

\* A tiny projected occurrence is not allowed to consume a mesh admission.
TinyNeverConsumesMeshBudget == mode[Tiny] \in {"box", "aggregate"}

\* Every candidate has at least structural coverage throughout this model.
CoverageNeverDisappears == \A candidate \in Candidates: mode[candidate] \in Modes

\* A terminal shortage is explicit and cannot conceal an affordable upgrade.
TerminalDebtExplicit ==
    /\ ~allocationPending /\ ~AllClassified
    => constrained

ConstrainedDebtHasNoAffordableUpgrade == constrained => ~AnyAffordable

SafeHeroEventuallyDirect ==
    [](profile = "safe" /\ allocationPending /\ mode["hero"] = "box" /\
       budget >= Cost("hero") => <>(mode["hero"] = "direct"))

LargeHeroEventuallyPop ==
    [](profile = "large" /\ allocationPending /\ mode["hero"] = "box" /\
       budget >= Cost("hero") => <>(mode["hero"] = "pop"))

EventuallyTerminal == <>(~allocationPending)

=============================================================================
