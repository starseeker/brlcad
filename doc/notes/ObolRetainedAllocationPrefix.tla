-------------------- MODULE ObolRetainedAllocationPrefix --------------------
\* Refinement contract between retained allocation and capacity search.
\*
\* One immutable scene revision supplies a deterministic visual-importance
\* order and a positive cost for every incremental upgrade.  Every budget
\* starts from the same coherent baseline and accepts upgrades in that order.
\* The first unaffordable upgrade terminates the allocation; later upgrades
\* may not fill the gap.  Consequently a larger budget can only extend the
\* population selected by a smaller budget.  Equal selected cost identifies
\* the same population, so capacity search may reuse a classified result
\* without presenting another frame.
\*
\* This model intentionally fixes the scene revision, baseline, importance
\* order, and upgrade costs.  Changing any of them is a new allocation
\* problem and must change the production transaction key.

EXTENDS Naturals, TLC

CONSTANT UpgradeCount, MaximumUpgradeCost, UpgradeCosts

ASSUME /\ UpgradeCount > 0
       /\ MaximumUpgradeCost > 0

Upgrades == 1..UpgradeCount
MaximumBudget == UpgradeCount * MaximumUpgradeCost

ASSUME UpgradeCosts \in [Upgrades -> 1..MaximumUpgradeCost]

\* A TLC fixture with an expensive first choice followed by cheap choices
\* exercises the gap which a skip-and-fill allocator handles incorrectly.
VariedUpgradeCosts ==
    [i \in Upgrades |->
        IF i = 1 THEN MaximumUpgradeCost
        ELSE 1 + ((i - 2) % MaximumUpgradeCost)]

RECURSIVE PrefixCost(_)
PrefixCost(count) ==
    IF count = 0
    THEN 0
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
    /\ lowerBudget \in 0..MaximumBudget
    /\ upperBudget \in lowerBudget..MaximumBudget
    /\ lowerPhase \in Phases
    /\ upperPhase \in Phases
    /\ lowerSelected \in 0..UpgradeCount
    /\ upperSelected \in 0..UpgradeCount
    /\ lowerSpent \in 0..MaximumBudget
    /\ upperSpent \in 0..MaximumBudget

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

NestedTerminalPopulations ==
    (lowerPhase = "done" /\ upperPhase = "done") =>
        lowerSelected <= upperSelected

EqualCostIdentifiesPopulation ==
    (lowerPhase = "done" /\ upperPhase = "done" /\
     lowerSpent = upperSpent) =>
        lowerSelected = upperSelected

Init ==
    /\ lowerBudget \in 0..MaximumBudget
    /\ upperBudget \in lowerBudget..MaximumBudget
    /\ lowerPhase = "selecting"
    /\ upperPhase = "selecting"
    /\ lowerSelected = 0
    /\ upperSelected = 0
    /\ lowerSpent = 0
    /\ upperSpent = 0

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
