-------------------------- MODULE ObolLodArbitration -------------------------
\* Bounded control-plane model for retained-allocation arbitration.
\*
\* This intentionally models only the policy contract, not rendering or
\* geometry.  "prominent" is a visible occurrence below its protected quality
\* floor.  The ordinary candidates represent otherwise useful refinements.
\* A transition consumes its coherent resident population's budget atomically.
\*
\* The production counterpart is MarginalUpgradeLess in
\* src/libBObol/retained_allocation_private.cpp.  It must not let optional
\* improvements overtake an affordable protected quality-floor transition;
\* when a budget cannot afford the next transition, the resulting debt is an
\* explicit, stable constrained outcome rather than an allocation loop.

EXTENDS Naturals, TLC

CONSTANT MaxEpoch

Candidates == {"prominent", "ordinary-a", "ordinary-b"}
Prominent == "prominent"
Ordinary == Candidates \ {Prominent}
InitialBudgets == 1..3

Cost(candidate) == IF candidate = Prominent THEN 2 ELSE 1

VARIABLES epoch,
          inputOpen,
          cuts,
          budget,
          allocationPending,
          constrained,
          previousEpoch,
          previousCuts

vars == <<epoch, inputOpen, cuts, budget, allocationPending, constrained,
          previousEpoch, previousCuts>>

TypeOK ==
    /\ epoch \in 0..MaxEpoch
    /\ inputOpen \in BOOLEAN
    /\ cuts \in [Candidates -> 0..1]
    /\ budget \in 0..3
    /\ allocationPending \in BOOLEAN
    /\ constrained \in BOOLEAN
    /\ previousEpoch \in 0..MaxEpoch
    /\ previousCuts \in [Candidates -> 0..1]

RememberPrevious ==
    /\ previousEpoch' = epoch
    /\ previousCuts' = cuts

Init ==
    /\ epoch = 0
    /\ inputOpen = FALSE
    /\ cuts = [candidate \in Candidates |-> 0]
    /\ budget \in InitialBudgets
    /\ allocationPending = TRUE
    /\ constrained = FALSE
    /\ previousEpoch = 0
    /\ previousCuts = cuts

\* A camera or source change invalidates allocation conclusions but not a
\* resident prefix.  In this small model every new epoch begins with no cuts,
\* so the monotonicity assertion deliberately applies only within one epoch.
BeginInput ==
    /\ ~inputOpen
    /\ ~allocationPending
    /\ epoch < MaxEpoch
    /\ epoch' = epoch + 1
    /\ inputOpen' = TRUE
    /\ cuts' = [candidate \in Candidates |-> 0]
    /\ budget' \in InitialBudgets
    /\ allocationPending' = FALSE
    /\ constrained' = FALSE
    /\ RememberPrevious

EndInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ allocationPending' = TRUE
    /\ UNCHANGED <<epoch, cuts, budget, constrained>>
    /\ RememberPrevious

AdmitProminent ==
    /\ allocationPending
    /\ ~inputOpen
    /\ cuts[Prominent] = 0
    /\ budget >= Cost(Prominent)
    /\ cuts' = [cuts EXCEPT ![Prominent] = 1]
    /\ budget' = budget - Cost(Prominent)
    /\ UNCHANGED <<epoch, inputOpen, allocationPending, constrained>>
    /\ RememberPrevious

\* Optional quality may be admitted only after the protected transition, or
\* after the allocator has proved that the protected transition is not
\* affordable in this epoch.
AdmitOrdinary(candidate) ==
    /\ candidate \in Ordinary
    /\ allocationPending
    /\ ~inputOpen
    /\ cuts[candidate] = 0
    /\ budget >= Cost(candidate)
    /\ cuts[Prominent] = 1 \/ budget < Cost(Prominent)
    /\ cuts' = [cuts EXCEPT ![candidate] = 1]
    /\ budget' = budget - Cost(candidate)
    /\ UNCHANGED <<epoch, inputOpen, allocationPending, constrained>>
    /\ RememberPrevious

AllCutsAdmitted == \A candidate \in Candidates: cuts[candidate] = 1
AnyAffordableCut ==
    \E candidate \in Candidates:
        /\ cuts[candidate] = 0
        /\ budget >= Cost(candidate)
        /\ (candidate = Prominent \/ cuts[Prominent] = 1 \/
            budget < Cost(Prominent))

Complete ==
    /\ allocationPending
    /\ AllCutsAdmitted
    /\ allocationPending' = FALSE
    /\ UNCHANGED <<epoch, inputOpen, cuts, budget, constrained>>
    /\ RememberPrevious

Constrain ==
    /\ allocationPending
    /\ ~AnyAffordableCut
    /\ ~AllCutsAdmitted
    /\ allocationPending' = FALSE
    /\ constrained' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, cuts, budget>>
    /\ RememberPrevious

Next ==
    \/ BeginInput
    \/ EndInput
    \/ AdmitProminent
    \/ \E candidate \in Ordinary: AdmitOrdinary(candidate)
    \/ Complete
    \/ Constrain

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(EndInput)
    /\ WF_vars(AdmitProminent)
    /\ \A candidate \in Ordinary: WF_vars(AdmitOrdinary(candidate))
    /\ WF_vars(Complete)
    /\ WF_vars(Constrain)

BudgetSafety == budget >= 0

\* An affordable protected candidate has not been displaced by an optional
\* candidate in the same allocation pass.
CoveragePriority ==
    cuts[Prominent] = 0 /\ budget >= Cost(Prominent)
    => \A candidate \in Ordinary: cuts[candidate] = 0

\* Within an unchanged view/source epoch, transitions only add coherent cuts.
SameEpochCutsMonotone ==
    previousEpoch = epoch =>
        \A candidate \in Candidates: previousCuts[candidate] <= cuts[candidate]

\* At a terminal shortfall the reason is explicit, and no policy-approved
\* transition remains affordable.  This prevents the balancing/refining cycle
\* from being encoded as an implicit terminal state.
TerminalDebtExplicit ==
    /\ ~inputOpen /\ ~allocationPending /\ ~AllCutsAdmitted
    => constrained

ConstrainedDebtHasNoAffordableUpgrade ==
    constrained => ~AnyAffordableCut

ProminentFloorNonStarvation ==
    [](allocationPending /\ ~inputOpen /\ cuts[Prominent] = 0 /\
       budget >= Cost(Prominent) => <>(cuts[Prominent] = 1))

EventuallyTerminal ==
    <>(~inputOpen /\ ~allocationPending)

=============================================================================
