------------------------- MODULE ObolResidentGrowth -------------------------
\* Availability discovered after a retained scene allocation is a typed
\* transaction.  A suffix requested by an active capacity candidate remains
\* inside that candidate; creating a second growth owner would make each wait
\* for the other.  Every independent coalesced growth edge first runs an
\* ordinary immutable-suffix drain and then enables exactly one scene-wide
\* allocation.  A growth edge arriving during that drain requires one more
\* drain.  Coverage is an output of the drain, not its admission guard: an
\* incomplete population consumes the growth transaction and transfers to one
\* ordinary coverage successor instead of restarting the same pre-drain pass.
\* Drain completion may never fall through to ordinary capacity calibration or
\* restart the just-completed pass while growth still owns the cursor.

EXTENDS Naturals, TLC

CONSTANT GrowthCount

ASSUME GrowthCount > 0

Phases == {"idle", "drain_required", "drain_active",
           "drain_active_dirty", "allocation_ready", "terminal"}

CoverageStates == {"complete", "incomplete"}

VARIABLES phase, growthRemaining, capacityOwnedGrowth, drainsStarted,
          drainsCompleted, allocations, coverageTransfers, coverage, cursor,
          capacityAvailable, capacityAbsorbedGrowth, lastAction

vars == <<phase, growthRemaining, capacityOwnedGrowth, drainsStarted,
          drainsCompleted, allocations, coverageTransfers, coverage, cursor,
          capacityAvailable, capacityAbsorbedGrowth, lastAction>>

TypeOK ==
    /\ phase \in Phases
    /\ growthRemaining \in 0..GrowthCount
    /\ capacityOwnedGrowth \in 0..GrowthCount
    /\ drainsStarted \in 0..GrowthCount
    /\ drainsCompleted \in 0..GrowthCount
    /\ allocations \in 0..GrowthCount
    /\ coverageTransfers \in 0..GrowthCount
    /\ coverage \in CoverageStates
    /\ cursor \in {"idle", "ordinary", "drain", "capacity"}
    /\ capacityAvailable \in BOOLEAN
    /\ capacityAbsorbedGrowth \in BOOLEAN
    /\ lastAction \in {"init", "growth", "capacity_growth",
                        "capacity_begin", "capacity_complete",
                        "capacity_restart", "yield",
                        "begin", "complete", "allocate",
                        "coverage_transfer", "coverage_complete", "finish"}

OrderedWork ==
    /\ allocations + coverageTransfers <= drainsCompleted
    /\ drainsCompleted <= drainsStarted

DrainCompletionHasTypedSuccessor ==
    lastAction = "complete" => phase \in {"allocation_ready", "drain_required"}

AllocationRequiresCompletedDrain ==
    lastAction = "allocate" => drainsCompleted > 0

CoverageTransferRequiresCompletedDrain ==
    lastAction = "coverage_transfer" => drainsCompleted > 0

CapacityRestartExcludesGrowth ==
    lastAction = "capacity_restart" => phase = "idle"

CapacityGrowthRetainsCandidate ==
    lastAction = "capacity_growth" =>
        /\ phase = "idle"
        /\ cursor = "capacity"

CapacityBeginCreatesProducer ==
    lastAction = "capacity_begin" =>
        /\ phase = "idle"
        /\ cursor = "capacity"
        /\ ~capacityAvailable

CapacityAbsorptionIsRecorded ==
    capacityAbsorbedGrowth => ~capacityAvailable

DrainCursorConsistent ==
    (cursor = "drain") <=>
        phase \in {"drain_active", "drain_active_dirty"}

Init ==
    /\ phase = "idle"
    /\ growthRemaining = GrowthCount
    /\ capacityOwnedGrowth = 0
    /\ drainsStarted = 0
    /\ drainsCompleted = 0
    /\ allocations = 0
    /\ coverageTransfers = 0
    /\ coverage \in CoverageStates
    \* Exercise both the independent-growth counterexample and a capacity
    \* candidate which directly owns the suffix it requested.
    /\ cursor \in {"ordinary", "capacity"}
    /\ capacityAvailable = (cursor = "ordinary")
    /\ capacityAbsorbedGrowth = FALSE
    /\ lastAction = "init"

PublishGrowth ==
    /\ growthRemaining > 0
    /\ growthRemaining' = growthRemaining - 1
    /\ IF cursor = "capacity"
          THEN /\ phase' = phase
               /\ capacityOwnedGrowth' = capacityOwnedGrowth + 1
               /\ lastAction' = "capacity_growth"
          ELSE /\ phase' =
                       CASE phase \in {"idle", "allocation_ready"} ->
                                "drain_required"
                          [] phase = "drain_active" -> "drain_active_dirty"
                          [] OTHER -> phase
               /\ capacityOwnedGrowth' = capacityOwnedGrowth
               /\ lastAction' = "growth"
    /\ UNCHANGED <<drainsStarted, drainsCompleted, allocations,
                    coverageTransfers, coverage, cursor,
                    capacityAvailable, capacityAbsorbedGrowth>>

BeginCapacity ==
    /\ capacityAvailable
    /\ cursor = "idle"
    /\ phase \in {"idle", "drain_required", "allocation_ready"}
    /\ phase' = "idle"
    /\ cursor' = "capacity"
    /\ capacityAvailable' = FALSE
    /\ capacityAbsorbedGrowth' = (phase # "idle")
    /\ lastAction' = "capacity_begin"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsStarted,
                    drainsCompleted, allocations, coverageTransfers, coverage>>

CompleteCapacity ==
    /\ cursor = "capacity"
    /\ phase = "idle"
    /\ cursor' = "idle"
    /\ lastAction' = "capacity_complete"
    /\ UNCHANGED <<phase, growthRemaining, capacityOwnedGrowth,
                    drainsStarted, drainsCompleted, allocations,
                    coverageTransfers, coverage,
                    capacityAvailable, capacityAbsorbedGrowth>>

CompleteOrdinary ==
    /\ cursor = "ordinary"
    /\ cursor' = IF phase = "idle" THEN "ordinary" ELSE "idle"
    /\ lastAction' = IF phase = "idle" THEN "capacity_restart" ELSE "yield"
    /\ UNCHANGED <<phase, growthRemaining, capacityOwnedGrowth,
                    drainsStarted, drainsCompleted, allocations,
                    coverageTransfers, coverage,
                    capacityAvailable, capacityAbsorbedGrowth>>

BeginDrain ==
    /\ phase = "drain_required"
    /\ cursor = "idle"
    /\ phase' = "drain_active"
    /\ cursor' = "drain"
    /\ drainsStarted' = drainsStarted + 1
    /\ lastAction' = "begin"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsCompleted, allocations,
                    coverageTransfers, coverage, capacityAvailable,
                    capacityAbsorbedGrowth>>

CompleteDrain ==
    /\ phase \in {"drain_active", "drain_active_dirty"}
    /\ cursor = "drain"
    /\ phase' = IF phase = "drain_active_dirty"
                   THEN "drain_required"
                   ELSE "allocation_ready"
    /\ cursor' = "idle"
    /\ drainsCompleted' = drainsCompleted + 1
    /\ lastAction' = "complete"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsStarted, allocations,
                    coverageTransfers, coverage, capacityAvailable,
                    capacityAbsorbedGrowth>>

Allocate ==
    /\ phase = "allocation_ready"
    /\ cursor = "idle"
    /\ coverage = "complete"
    /\ phase' = "idle"
    /\ allocations' = allocations + 1
    /\ lastAction' = "allocate"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsStarted, drainsCompleted,
                    coverageTransfers, coverage, cursor, capacityAvailable,
                    capacityAbsorbedGrowth>>

TransferIncompleteCoverage ==
    /\ phase = "allocation_ready"
    /\ cursor = "idle"
    /\ coverage = "incomplete"
    /\ phase' = "idle"
    /\ cursor' = "ordinary"
    /\ coverageTransfers' = coverageTransfers + 1
    /\ lastAction' = "coverage_transfer"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsStarted, drainsCompleted,
                    allocations, coverage, capacityAvailable,
                    capacityAbsorbedGrowth>>

CompleteTransferredCoverage ==
    /\ phase = "idle"
    /\ cursor = "ordinary"
    /\ coverage = "incomplete"
    /\ cursor' = "idle"
    /\ coverage' = "complete"
    /\ lastAction' = "coverage_complete"
    /\ UNCHANGED <<phase, growthRemaining, capacityOwnedGrowth, drainsStarted, drainsCompleted,
                    allocations, coverageTransfers, capacityAvailable,
                    capacityAbsorbedGrowth>>

Finish ==
    /\ phase = "idle"
    /\ growthRemaining = 0
    /\ cursor = "idle"
    /\ ~capacityAvailable
    /\ phase' = "terminal"
    /\ lastAction' = "finish"
    /\ UNCHANGED <<growthRemaining, capacityOwnedGrowth, drainsStarted, drainsCompleted,
                    allocations, coverageTransfers, coverage, cursor,
                    capacityAvailable, capacityAbsorbedGrowth>>

Next == PublishGrowth \/ BeginCapacity \/ CompleteCapacity \/ CompleteOrdinary \/
        BeginDrain \/ CompleteDrain \/
        Allocate \/ TransferIncompleteCoverage \/
        CompleteTransferredCoverage \/ Finish

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishGrowth)
    /\ WF_vars(BeginCapacity)
    /\ WF_vars(CompleteCapacity)
    /\ WF_vars(CompleteOrdinary)
    /\ WF_vars(BeginDrain)
    /\ WF_vars(CompleteDrain)
    /\ WF_vars(Allocate)
    /\ WF_vars(TransferIncompleteCoverage)
    /\ WF_vars(CompleteTransferredCoverage)
    /\ WF_vars(Finish)

EventuallySettles == <> (phase = "terminal")

=============================================================================
