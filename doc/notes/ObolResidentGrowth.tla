------------------------- MODULE ObolResidentGrowth -------------------------
\* Availability discovered after a retained scene allocation is a typed
\* transaction.  Every coalesced growth edge must first run an ordinary
\* immutable-suffix drain and then enable exactly one scene-wide allocation.
\* A growth edge arriving during that drain requires one more drain before the
\* allocation.  Drain completion may never fall through to ordinary capacity
\* calibration or restart the just-completed pass.

EXTENDS Naturals, TLC

CONSTANT GrowthCount

ASSUME GrowthCount > 0

Phases == {"idle", "drain_required", "drain_active",
           "drain_active_dirty", "allocation_ready", "terminal"}

VARIABLES phase, growthRemaining, drainsStarted, drainsCompleted,
          allocations, cursor, lastAction

vars == <<phase, growthRemaining, drainsStarted, drainsCompleted,
          allocations, cursor, lastAction>>

TypeOK ==
    /\ phase \in Phases
    /\ growthRemaining \in 0..GrowthCount
    /\ drainsStarted \in 0..GrowthCount
    /\ drainsCompleted \in 0..GrowthCount
    /\ allocations \in 0..GrowthCount
    /\ cursor \in {"idle", "ordinary", "drain"}
    /\ lastAction \in {"init", "growth", "capacity_restart", "yield",
                        "begin", "complete", "allocate", "finish"}

OrderedWork == allocations <= drainsCompleted /\ drainsCompleted <= drainsStarted

DrainCompletionHasTypedSuccessor ==
    lastAction = "complete" => phase \in {"allocation_ready", "drain_required"}

AllocationRequiresCompletedDrain ==
    lastAction = "allocate" => drainsCompleted > 0

CapacityRestartExcludesGrowth ==
    lastAction = "capacity_restart" => phase = "idle"

DrainCursorConsistent ==
    (cursor = "drain") <=>
        phase \in {"drain_active", "drain_active_dirty"}

Init ==
    /\ phase = "idle"
    /\ growthRemaining = GrowthCount
    /\ drainsStarted = 0
    /\ drainsCompleted = 0
    /\ allocations = 0
    \* Model a capacity-blocked ordinary pass already in flight when an
    \* immutable suffix becomes resident.  This is the runtime counterexample
    \* which used to restart forever instead of yielding its cursor.
    /\ cursor = "ordinary"
    /\ lastAction = "init"

PublishGrowth ==
    /\ growthRemaining > 0
    /\ growthRemaining' = growthRemaining - 1
    /\ phase' = CASE phase \in {"idle", "allocation_ready"} ->
                        "drain_required"
                   [] phase = "drain_active" -> "drain_active_dirty"
                   [] OTHER -> phase
    /\ lastAction' = "growth"
    /\ UNCHANGED <<drainsStarted, drainsCompleted, allocations, cursor>>

CompleteOrdinary ==
    /\ cursor = "ordinary"
    /\ cursor' = IF phase = "idle" THEN "ordinary" ELSE "idle"
    /\ lastAction' = IF phase = "idle" THEN "capacity_restart" ELSE "yield"
    /\ UNCHANGED <<phase, growthRemaining, drainsStarted, drainsCompleted,
                    allocations>>

BeginDrain ==
    /\ phase = "drain_required"
    /\ cursor = "idle"
    /\ phase' = "drain_active"
    /\ cursor' = "drain"
    /\ drainsStarted' = drainsStarted + 1
    /\ lastAction' = "begin"
    /\ UNCHANGED <<growthRemaining, drainsCompleted, allocations>>

CompleteDrain ==
    /\ phase \in {"drain_active", "drain_active_dirty"}
    /\ cursor = "drain"
    /\ phase' = IF phase = "drain_active_dirty"
                   THEN "drain_required"
                   ELSE "allocation_ready"
    /\ cursor' = "idle"
    /\ drainsCompleted' = drainsCompleted + 1
    /\ lastAction' = "complete"
    /\ UNCHANGED <<growthRemaining, drainsStarted, allocations>>

Allocate ==
    /\ phase = "allocation_ready"
    /\ cursor = "idle"
    /\ phase' = "idle"
    /\ allocations' = allocations + 1
    /\ lastAction' = "allocate"
    /\ UNCHANGED <<growthRemaining, drainsStarted, drainsCompleted, cursor>>

Finish ==
    /\ phase = "idle"
    /\ growthRemaining = 0
    /\ cursor = "idle"
    /\ phase' = "terminal"
    /\ lastAction' = "finish"
    /\ UNCHANGED <<growthRemaining, drainsStarted, drainsCompleted,
                    allocations, cursor>>

Next == PublishGrowth \/ CompleteOrdinary \/ BeginDrain \/ CompleteDrain \/
        Allocate \/ Finish

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishGrowth)
    /\ WF_vars(CompleteOrdinary)
    /\ WF_vars(BeginDrain)
    /\ WF_vars(CompleteDrain)
    /\ WF_vars(Allocate)
    /\ WF_vars(Finish)

EventuallySettles == <> (phase = "terminal")

=============================================================================
