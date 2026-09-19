--------------------------- MODULE ObolLodAdmission -------------------------
\* Bounded representation-admission contract for one immutable view demand.
\*
\* Discovery and mesh availability are intentionally concurrent.  A cold
\* source may therefore publish a useful PoP mesh before the complete source
\* inventory proves that the scene is small enough for direct terminal
\* geometry.  Certification must be allowed to promote that existing PoP
\* representation atomically, charging only the replacement's marginal cost.
\*
\* Structural boxes are truthful startup coverage, never a terminal
\* representation.  A significant occurrence that cannot afford even its
\* minimum mesh may become a persistent proxy only with an explicit capacity
\* witness.  Tiny occurrences use the aggregate channel and consume no
\* per-occurrence mesh budget.

EXTENDS Naturals, TLC

\* Source profile certification, not a workload/cardinality profile.
Profiles == {"unknown", "safe", "large"}
Candidates == {"hero", "ordinary", "tiny"}
Significant == Candidates \ {"tiny"}
Modes == {"box", "aggregate", "proxy", "pop", "direct"}
InitialCapacities == 0..4

PopCost(candidate) == IF candidate = "hero" THEN 2 ELSE 1
DirectCost(candidate) == 2

ModeCost(candidate, candidateMode) ==
    IF candidateMode = "pop" THEN PopCost(candidate)
    ELSE IF candidateMode = "direct" THEN DirectCost(candidate)
    ELSE 0

TotalCost(candidateModes) ==
    ModeCost("hero", candidateModes["hero"]) +
    ModeCost("ordinary", candidateModes["ordinary"])

Replace(candidateModes, candidate, replacement) ==
    [candidateModes EXCEPT ![candidate] = replacement]

VARIABLES profile,
          capacity,
          mode,
          planning,
          constraintWitness

vars == <<profile, capacity, mode, planning, constraintWitness>>

TypeOK ==
    /\ profile \in Profiles
    /\ capacity \in InitialCapacities
    /\ mode \in [Candidates -> Modes]
    /\ planning \in BOOLEAN
    /\ constraintWitness \in BOOLEAN

Init ==
    /\ profile = "unknown"
    /\ capacity \in InitialCapacities
    /\ mode = [candidate \in Candidates |-> "box"]
    /\ planning = TRUE
    /\ constraintWitness = FALSE

ClassifyTiny ==
    /\ planning
    /\ mode["tiny"] = "box"
    /\ mode' = [mode EXCEPT !["tiny"] = "aggregate"]
    /\ UNCHANGED <<profile, capacity, planning, constraintWitness>>

Discover(nextProfile) ==
    /\ planning
    /\ profile = "unknown"
    /\ nextProfile \in {"safe", "large"}
    /\ profile' = nextProfile
    /\ UNCHANGED <<capacity, mode, planning, constraintWitness>>

AdmitPop(candidate) ==
    LET replacement == Replace(mode, candidate, "pop") IN
    /\ planning
    /\ candidate \in Significant
    /\ mode[candidate] = "box"
    /\ TotalCost(replacement) <= capacity
    /\ mode' = replacement
    /\ UNCHANGED <<profile, capacity, planning, constraintWitness>>

AdmitDirect(candidate) ==
    LET replacement == Replace(mode, candidate, "direct") IN
    /\ planning
    /\ profile = "safe"
    /\ candidate \in Significant
    /\ mode[candidate] \in {"box", "pop"}
    /\ TotalCost(replacement) <= capacity
    /\ mode' = replacement
    /\ UNCHANGED <<profile, capacity, planning, constraintWitness>>

PublishPersistentProxy(candidate) ==
    LET popReplacement == Replace(mode, candidate, "pop")
        directReplacement == Replace(mode, candidate, "direct") IN
    /\ planning
    /\ profile # "unknown"
    /\ candidate \in Significant
    /\ mode[candidate] = "box"
    /\ TotalCost(popReplacement) > capacity
    /\ (profile # "safe" \/ TotalCost(directReplacement) > capacity)
    /\ mode' = Replace(mode, candidate, "proxy")
    /\ constraintWitness' = TRUE
    /\ UNCHANGED <<profile, capacity, planning>>

AllPresented == \A candidate \in Candidates: mode[candidate] # "box"

AffordableDesiredUpgrade(candidate) ==
    IF profile = "safe"
    THEN /\ mode[candidate] \in {"box", "pop"}
         /\ TotalCost(Replace(mode, candidate, "direct")) <= capacity
    ELSE /\ mode[candidate] = "box"
         /\ TotalCost(Replace(mode, candidate, "pop")) <= capacity

AnyAffordableDesiredUpgrade ==
    \E candidate \in Significant: AffordableDesiredUpgrade(candidate)

Complete ==
    /\ planning
    /\ profile # "unknown"
    /\ AllPresented
    /\ ~AnyAffordableDesiredUpgrade
    /\ planning' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, constraintWitness>>

Next ==
    \/ ClassifyTiny
    \/ \E nextProfile \in {"safe", "large"}: Discover(nextProfile)
    \/ \E candidate \in Significant: AdmitPop(candidate)
    \/ \E candidate \in Significant: AdmitDirect(candidate)
    \/ \E candidate \in Significant: PublishPersistentProxy(candidate)
    \/ Complete

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ClassifyTiny)
    /\ \A nextProfile \in {"safe", "large"}:
           WF_vars(Discover(nextProfile))
    /\ \A candidate \in Significant: WF_vars(AdmitPop(candidate))
    /\ \A candidate \in Significant: WF_vars(AdmitDirect(candidate))
    /\ \A candidate \in Significant:
           WF_vars(PublishPersistentProxy(candidate))
    /\ WF_vars(Complete)

BudgetSafety == TotalCost(mode) <= capacity

TinyNeverConsumesMeshBudget == mode["tiny"] \in {"box", "aggregate"}

\* A structural startup representation cannot be mistaken for convergence.
TerminalHasNoStructuralBoxes ==
    ~planning => \A candidate \in Candidates: mode[candidate] # "box"

\* A persistent proxy is a certified constrained choice, never an admission
\* shortcut taken merely because discovery or source preparation is pending.
PersistentProxyHasWitness ==
    (\E candidate \in Significant: mode[candidate] = "proxy") =>
        constraintWitness

\* Once published, a richer representation cannot regress during the same
\* immutable demand.  New view/policy/capacity revisions are outside this
\* focused model and begin a new plan.
DirectIsAbsorbing ==
    \A candidate \in Significant:
        [](mode[candidate] = "direct" => [](
            mode[candidate] = "direct"))

\* This is the cold-ordering regression boundary.  Once complete discovery
\* certifies a safe scene, an affordable exact replacement remains reachable
\* even when the occurrence already owns a useful PoP mesh.
ZeroMarginalSafePopEventuallyDirect ==
    \A candidate \in Significant:
        [](profile = "safe" /\ planning /\ mode[candidate] = "pop" /\
           TotalCost(Replace(mode, candidate, "direct")) = TotalCost(mode)
           => <>(mode[candidate] = "direct"))

SafeFullCapacityEventuallyDirect ==
    [](profile = "safe" /\ planning /\ capacity = 4
       => <>(mode["hero"] = "direct" /\
              mode["ordinary"] = "direct"))

EventuallyTerminal == <>(~planning)

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => ~planning

=============================================================================
