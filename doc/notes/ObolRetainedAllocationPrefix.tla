-------------------- MODULE ObolRetainedAllocationPrefix --------------------
\* Refinement contract between retained allocation and capacity search.
\*
\* One immutable scene revision supplies a deterministic visual-importance
\* order and a non-negative cost for every incremental upgrade.  Every budget
\* starts from the same coherent baseline and accepts upgrades in that order.
\* The baseline is the exact presentation selected by the protected floor:
\* fixed records plus either the protected mesh cut or the aggregate point
\* proxy for every candidate.  A proxy replaces its mesh cost; it does not
\* erase its own presentation cost from the certificate.
\* The first unaffordable upgrade terminates the allocation; later upgrades
\* may not fill the gap.  Consequently a larger budget can only extend the
\* population selected by a smaller budget.  Zero-cost precision upgrades are
\* always consumed before a population is terminal.  Equal selected cost thus
\* identifies the same canonical population, so capacity search may reuse a
\* classified result without presenting another frame.
\*
\* This model intentionally fixes the scene revision, baseline, importance
\* order, and upgrade costs.  Changing any of them is a new allocation
\* problem and must change the production transaction key.

EXTENDS Naturals, TLC

CONSTANT UpgradeCount, MaximumUpgradeCost, UpgradeCosts, FixedCost,
         ProtectedMeshCosts, PointProxyCosts, PointEligible

ASSUME /\ UpgradeCount > 0
       /\ MaximumUpgradeCost > 0
       /\ FixedCost \in 0..MaximumUpgradeCost

Upgrades == 1..UpgradeCount

ASSUME UpgradeCosts \in [Upgrades -> 0..MaximumUpgradeCost]
ASSUME ProtectedMeshCosts \in [Upgrades -> 0..MaximumUpgradeCost]
ASSUME PointProxyCosts \in [Upgrades -> 0..MaximumUpgradeCost]
ASSUME PointEligible \in [Upgrades -> BOOLEAN]

\* A TLC fixture with an expensive first choice followed by cheap choices
\* exercises the gap which a skip-and-fill allocator handles incorrectly.
VariedUpgradeCosts ==
    [i \in Upgrades |->
        IF i = 1 THEN MaximumUpgradeCost
        ELSE IF i % 3 = 0 THEN 0
        ELSE 1 + ((i - 2) % MaximumUpgradeCost)]

VariedProtectedMeshCosts ==
    [i \in Upgrades |-> 1 + ((i - 1) % MaximumUpgradeCost)]

VariedPointProxyCosts ==
    [i \in Upgrades |-> IF i % 2 = 0 THEN 1 ELSE 0]

VariedPointEligibility ==
    [i \in Upgrades |-> i % 2 = 0]

RECURSIVE BaselinePrefixCost(_)
BaselinePrefixCost(count) ==
    IF count = 0
    THEN FixedCost
    ELSE BaselinePrefixCost(count - 1) +
        IF PointEligible[count]
        THEN PointProxyCosts[count]
        ELSE ProtectedMeshCosts[count]

BaselineCost == BaselinePrefixCost(UpgradeCount)
MaximumBudget == BaselineCost + UpgradeCount * MaximumUpgradeCost

RECURSIVE PrefixCost(_)
PrefixCost(count) ==
    IF count = 0
    THEN BaselineCost
    ELSE PrefixCost(count - 1) + UpgradeCosts[count]

Phases == {"selecting", "done"}

VARIABLES lowerBudget,
          upperBudget,
          lowerPhase,
          upperPhase,
          lowerSelected,
          upperSelected,
          lowerSpent,
          upperSpent

vars == <<lowerBudget, upperBudget, lowerPhase, upperPhase,
          lowerSelected, upperSelected, lowerSpent, upperSpent>>

TypeOK ==
    /\ lowerBudget \in BaselineCost..MaximumBudget
    /\ upperBudget \in lowerBudget..MaximumBudget
    /\ lowerPhase \in Phases
    /\ upperPhase \in Phases
    /\ lowerSelected \in 0..UpgradeCount
    /\ upperSelected \in 0..UpgradeCount
    /\ lowerSpent \in BaselineCost..MaximumBudget
    /\ upperSpent \in BaselineCost..MaximumBudget

SpendMatchesPrefix ==
    /\ lowerSpent = PrefixCost(lowerSelected)
    /\ upperSpent = PrefixCost(upperSelected)

BudgetSafe ==
    /\ lowerSpent <= lowerBudget
    /\ upperSpent <= upperBudget

TerminalAtFirstGap ==
    /\ lowerPhase = "done" =>
          \/ lowerSelected = UpgradeCount
          \/ lowerSpent + UpgradeCosts[lowerSelected + 1] > lowerBudget
    /\ upperPhase = "done" =>
          \/ upperSelected = UpgradeCount
          \/ upperSpent + UpgradeCosts[upperSelected + 1] > upperBudget

TerminalIncludesFreeSuccessor ==
    /\ (lowerPhase = "done" /\ lowerSelected < UpgradeCount =>
            UpgradeCosts[lowerSelected + 1] > 0)
    /\ (upperPhase = "done" /\ upperSelected < UpgradeCount =>
            UpgradeCosts[upperSelected + 1] > 0)

NestedTerminalPopulations ==
    (lowerPhase = "done" /\ upperPhase = "done") =>
        lowerSelected <= upperSelected

EqualCostIdentifiesPopulation ==
    (lowerPhase = "done" /\ upperPhase = "done" /\
     lowerSpent = upperSpent) =>
        lowerSelected = upperSelected

Init ==
    /\ lowerBudget \in BaselineCost..MaximumBudget
    /\ upperBudget \in lowerBudget..MaximumBudget
    /\ lowerPhase = "selecting"
    /\ upperPhase = "selecting"
    /\ lowerSelected = 0
    /\ upperSelected = 0
    /\ lowerSpent = BaselineCost
    /\ upperSpent = BaselineCost

AdvanceLower ==
    /\ lowerPhase = "selecting"
    /\ IF lowerSelected = UpgradeCount
          THEN /\ lowerPhase' = "done"
               /\ UNCHANGED <<lowerSelected, lowerSpent>>
          ELSE IF lowerSpent + UpgradeCosts[lowerSelected + 1] <= lowerBudget
               THEN /\ lowerSelected' = lowerSelected + 1
                    /\ lowerSpent' = lowerSpent +
                           UpgradeCosts[lowerSelected + 1]
                    /\ UNCHANGED lowerPhase
               ELSE /\ lowerPhase' = "done"
                    /\ UNCHANGED <<lowerSelected, lowerSpent>>
    /\ UNCHANGED <<lowerBudget, upperBudget, upperPhase,
                    upperSelected, upperSpent>>

AdvanceUpper ==
    /\ upperPhase = "selecting"
    /\ IF upperSelected = UpgradeCount
          THEN /\ upperPhase' = "done"
               /\ UNCHANGED <<upperSelected, upperSpent>>
          ELSE IF upperSpent + UpgradeCosts[upperSelected + 1] <= upperBudget
               THEN /\ upperSelected' = upperSelected + 1
                    /\ upperSpent' = upperSpent +
                           UpgradeCosts[upperSelected + 1]
                    /\ UNCHANGED upperPhase
               ELSE /\ upperPhase' = "done"
                    /\ UNCHANGED <<upperSelected, upperSpent>>
    /\ UNCHANGED <<lowerBudget, upperBudget, lowerPhase,
                    lowerSelected, lowerSpent>>

Next == AdvanceLower \/ AdvanceUpper

Spec == Init /\ [][Next]_vars
        /\ WF_vars(AdvanceLower)
        /\ WF_vars(AdvanceUpper)

EventuallyBothDone == <>(lowerPhase = "done" /\ upperPhase = "done")

=============================================================================
