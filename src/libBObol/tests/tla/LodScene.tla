------------------------------- MODULE LodScene -------------------------------
EXTENDS LodCommon

CONSTANTS
    MaxEpoch,
    ResidentLimit,
    InitialRenderBudget,
    HardQualityBudget,
    MaxRenderBudget,
    InitialMode,
    InjectStaleAdmission,
    InjectPrematureStable

Objects == {"A", "B", "C"}
Cuts == 0..2
CutDomain == -1..2
CapacityStates == {"Headroom", "NearTarget", "Overloaded"}
InitialModes == {"Normal", "Pressure", "Cancelled", "Overloaded"}

ASSUME
    /\ MaxEpoch \in Nat \ {0}
    /\ ResidentLimit \in Nat \ {0}
    /\ InitialRenderBudget \in 3..MaxRenderBudget
    /\ HardQualityBudget \in 3..MaxRenderBudget
    /\ MaxRenderBudget \in Nat
    /\ InitialMode \in InitialModes
    /\ InjectStaleAdmission \in BOOLEAN
    /\ InjectPrematureStable \in BOOLEAN

VARIABLES
    quiescent,
    sourceActive,
    viewEpoch,
    policyEpoch,
    inventoryEpoch,
    routingId,
    generation,
    interactive,
    gestureActive,
    quietPending,
    covered,
    activeCut,
    residentCut,
    requestedCut,
    allocatedCut,
    allocationView,
    allocationPolicy,
    scanOrder,
    memoryLimited,
    tasks,
    results,
    publicationPending,
    framePending,
    renderBudget,
    capacity,
    calibrationPending,
    budgetRecoveryPending,
    cpuPressure,
    gpuPressure,
    recoveryPending,
    compactionPending,
    compactionPlanning,
    staleAdmissionCount

vars == <<
    quiescent, sourceActive, viewEpoch, policyEpoch, inventoryEpoch, routingId,
    generation, interactive, gestureActive, quietPending, covered, activeCut,
    residentCut, requestedCut, allocatedCut, allocationView, allocationPolicy,
    scanOrder, memoryLimited, tasks, results, publicationPending, framePending,
    renderBudget, capacity, calibrationPending, budgetRecoveryPending,
    cpuPressure, gpuPressure, recoveryPending, compactionPending,
    compactionPlanning, staleAdmissionCount
>>

identityState == <<
    sourceActive, viewEpoch, policyEpoch, inventoryEpoch, routingId, generation
>>
interactionState == <<interactive, gestureActive, quietPending>>
geometryState == <<
    covered, activeCut, residentCut, requestedCut, allocatedCut,
    allocationView, allocationPolicy, scanOrder, memoryLimited
>>
workState == <<tasks, results, publicationPending, framePending>>
budgetState == <<
    renderBudget, capacity, calibrationPending, budgetRecoveryPending
>>
pressureState == <<
    cpuPressure, gpuPressure, recoveryPending, compactionPending,
    compactionPlanning
>>

Orders == {
    <<"A", "B", "C">>, <<"A", "C", "B">>,
    <<"B", "A", "C">>, <<"B", "C", "A">>,
    <<"C", "A", "B">>, <<"C", "B", "A">>
}

Importance == [A |-> 3, B |-> 2, C |-> 1]
ProtectedCut == [A |-> 1, B |-> 0, C |-> 0]
InitialRequestedCut == [A |-> 2, B |-> 2, C |-> 2]

RenderAt(c) ==
    CASE c = -1 -> 0
      [] c = 0 -> 1
      [] c = 1 -> 3
      [] c = 2 -> 6

MemoryAt(c) ==
    CASE c = -1 -> 0
      [] c = 0 -> 2
      [] c = 1 -> 5
      [] c = 2 -> 9

QualityAt(c) ==
    CASE c = 0 -> 0
      [] c = 1 -> 10
      [] c = 2 -> 15

TotalRender(cuts) ==
    RenderAt(cuts["A"]) + RenderAt(cuts["B"]) + RenderAt(cuts["C"])

TotalMemory(cuts) ==
    MemoryAt(cuts["A"]) + MemoryAt(cuts["B"]) + MemoryAt(cuts["C"])

AllocationScore(cuts) ==
    Importance["A"] * QualityAt(cuts["A"])
    + Importance["B"] * QualityAt(cuts["B"])
    + Importance["C"] * QualityAt(cuts["C"])

CutMaps == [Objects -> Cuts]

ProtectedFloorCost == TotalRender(ProtectedCut)
ProtectedFloorFeasible == ProtectedFloorCost <= HardQualityBudget

EffectiveBudgetFor(budget) ==
    IF ProtectedFloorFeasible
    THEN MaxNat(budget, ProtectedFloorCost)
    ELSE budget

EffectiveRenderBudget == EffectiveBudgetFor(renderBudget)

FeasibleAllocation(cuts) ==
    /\ cuts \in CutMaps
    /\ \A o \in Objects: cuts[o] <= requestedCut[o]
    /\ TotalRender(cuts) <= EffectiveRenderBudget
    /\ (ProtectedFloorFeasible =>
            \A o \in Objects: cuts[o] >= ProtectedCut[o])

BestAllocations == {
    candidate \in CutMaps:
        FeasibleAllocation(candidate)
        /\ \A alternative \in CutMaps:
            FeasibleAllocation(alternative) =>
                AllocationScore(candidate) >= AllocationScore(alternative)
}

AllocationCurrent ==
    sourceActive
    /\ allocationView = viewEpoch
    /\ allocationPolicy = policyEpoch
    /\ allocationView # 0
    /\ allocationPolicy # 0

NoRequest == [
    valid |-> FALSE,
    generation |-> 0,
    viewEpoch |-> 0,
    policyEpoch |-> 0,
    inventoryEpoch |-> 0,
    routingId |-> 0,
    cut |-> -1
]

MakeRequest(g, v, p, i, r, c) == [
    valid |-> TRUE,
    generation |-> g,
    viewEpoch |-> v,
    policyEpoch |-> p,
    inventoryEpoch |-> i,
    routingId |-> r,
    cut |-> c
]

RequestType(request) ==
    /\ request.valid \in BOOLEAN
    /\ request.generation \in 0..MaxEpoch
    /\ request.viewEpoch \in 0..MaxEpoch
    /\ request.policyEpoch \in 0..MaxEpoch
    /\ request.inventoryEpoch \in 0..MaxEpoch
    /\ request.routingId \in 0..MaxEpoch
    /\ request.cut \in CutDomain

RequestCurrent(request) ==
    request.valid
    /\ sourceActive
    /\ request.generation = generation
    /\ request.viewEpoch = viewEpoch
    /\ request.policyEpoch = policyEpoch
    /\ request.inventoryEpoch = inventoryEpoch
    /\ request.routingId = routingId

HasTask == \E o \in Objects: tasks[o].valid
HasResult == \E o \in Objects: results[o].valid

ResidentWith(o, cut) ==
    [residentCut EXCEPT ![o] = MaxNat(@, cut)]

CanAdmit(o, cut) ==
    cut = 0 \/ TotalMemory(ResidentWith(o, cut)) <= ResidentLimit

CurrentVisualSettled ==
    \A o \in Objects:
        activeCut[o] = allocatedCut[o] \/ o \in memoryLimited

ForegroundWitnesses ==
    (IF sourceActive /\ covered # Objects THEN {"CoverageCursor"} ELSE {})
    \cup (IF sourceActive /\ covered = Objects /\ ~AllocationCurrent
        THEN {"Allocation"} ELSE {})
    \cup (IF AllocationCurrent /\ ~CurrentVisualSettled
        THEN {"Application"} ELSE {})
    \cup (IF HasTask THEN {"ServiceTask"} ELSE {})
    \cup (IF HasResult THEN {"ServiceResult"} ELSE {})
    \cup (IF publicationPending THEN {"Publication"} ELSE {})
    \cup (IF framePending THEN {"PresentationFrame"} ELSE {})
    \cup (IF gestureActive THEN {"InteractionInput"} ELSE {})
    \cup (IF interactive \/ quietPending THEN {"QuietTimer"} ELSE {})
    \cup (IF calibrationPending THEN {"Calibration"} ELSE {})
    \cup (IF budgetRecoveryPending THEN {"BudgetRecovery"} ELSE {})
    \cup (IF recoveryPending THEN {"PressureRecovery"} ELSE {})
    \cup (IF compactionPlanning THEN {"CompactionCursor"} ELSE {})

CoordinatorPhase ==
    IF interactive \/ gestureActive THEN Interactive
    ELSE IF sourceActive /\ covered # Objects THEN Coverage
    ELSE IF ~sourceActive THEN Fallback
    ELSE IF compactionPlanning \/ recoveryPending THEN Compacting
    ELSE IF ForegroundWitnesses # {} THEN Settling
    ELSE Stable

\* The mutation branch is exercised only by an expected-failure test.  It
\* intentionally omits the completed-presentation barrier.
ReportedStable ==
    IF InjectPrematureStable
    THEN sourceActive
        /\ covered = Objects
        /\ AllocationCurrent
        /\ CurrentVisualSettled
        /\ ~interactive
        /\ ~gestureActive
        /\ ~quietPending
        /\ ~HasTask
        /\ ~HasResult
        /\ ~publicationPending
        /\ ~calibrationPending
        /\ ~budgetRecoveryPending
        /\ ~recoveryPending
        /\ ~compactionPlanning
    ELSE CoordinatorPhase = Stable
        /\ AllocationCurrent
        /\ CurrentVisualSettled

FallbackSettled ==
    ~sourceActive /\ ForegroundWitnesses = {} /\ ~compactionPending

Terminal == ReportedStable \/ FallbackSettled

PerformanceLimited ==
    ReportedStable /\ \E o \in Objects: allocatedCut[o] < requestedCut[o]

MemoryLimitedTerminal ==
    ReportedStable /\
        (memoryLimited # {} \/ cpuPressure \/ gpuPressure
         \/ TotalMemory(residentCut) > ResidentLimit)

FullRequestedQuality ==
    ReportedStable /\
        \A o \in Objects: activeCut[o] >= requestedCut[o]

Init ==
    /\ quiescent = FALSE
    /\ sourceActive = (InitialMode # "Cancelled")
    /\ viewEpoch = 1
    /\ policyEpoch = 1
    /\ inventoryEpoch = 1
    /\ routingId = 1
    /\ generation = IF InitialMode = "Cancelled" THEN 2 ELSE 1
    /\ interactive = FALSE
    /\ gestureActive = FALSE
    /\ quietPending = FALSE
    /\ covered = {}
    /\ activeCut = [o \in Objects |-> -1]
    /\ residentCut = [o \in Objects |-> -1]
    /\ requestedCut = InitialRequestedCut
    /\ allocatedCut = [o \in Objects |-> -1]
    /\ allocationView = 0
    /\ allocationPolicy = 0
    /\ scanOrder = <<"A", "B", "C">>
    /\ memoryLimited = {}
    /\ tasks = [o \in Objects |-> NoRequest]
    /\ results = [o \in Objects |->
            IF InitialMode = "Cancelled" /\ o = "A"
            THEN MakeRequest(1, 1, 1, 1, 1, 0)
            ELSE NoRequest]
    /\ publicationPending = FALSE
    /\ framePending = FALSE
    /\ renderBudget = InitialRenderBudget
    /\ capacity = IF InitialMode = "Overloaded"
            THEN "Overloaded" ELSE "NearTarget"
    /\ calibrationPending = (InitialMode = "Overloaded")
    /\ budgetRecoveryPending = FALSE
    /\ cpuPressure = (InitialMode = "Pressure")
    /\ gpuPressure = FALSE
    /\ recoveryPending = (InitialMode = "Pressure")
    /\ compactionPending = (InitialMode = "Pressure")
    /\ compactionPlanning = FALSE
    /\ staleAdmissionCount = 0

QuiesceEnvironment ==
    /\ ~quiescent
    /\ quiescent' = TRUE
    /\ UNCHANGED <<
        identityState, interactionState, geometryState, workState, budgetState,
        pressureState, staleAdmissionCount
       >>

InvalidateView ==
    /\ ~quiescent
    /\ sourceActive
    /\ viewEpoch < MaxEpoch
    /\ viewEpoch' = viewEpoch + 1
    /\ sourceActive' = sourceActive
    /\ policyEpoch' = policyEpoch
    /\ inventoryEpoch' = inventoryEpoch
    /\ routingId' = routingId
    /\ generation' = generation
    /\ covered' = {}
    /\ activeCut' = activeCut
    /\ residentCut' = residentCut
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = 0
    /\ allocationPolicy' = 0
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = {}
    /\ compactionPlanning' = FALSE
    /\ compactionPending' = compactionPending \/ recoveryPending
    /\ cpuPressure' = cpuPressure
    /\ gpuPressure' = gpuPressure
    /\ recoveryPending' = recoveryPending
    /\ UNCHANGED <<
        quiescent, interactionState, workState, budgetState,
        staleAdmissionCount
       >>

ChangePolicy ==
    /\ ~quiescent
    /\ sourceActive
    /\ policyEpoch < MaxEpoch
    /\ policyEpoch' = policyEpoch + 1
    /\ sourceActive' = sourceActive
    /\ viewEpoch' = viewEpoch
    /\ inventoryEpoch' = inventoryEpoch
    /\ routingId' = routingId
    /\ generation' = generation
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = 0
    /\ allocationPolicy' = 0
    /\ covered' = covered
    /\ activeCut' = activeCut
    /\ residentCut' = residentCut
    /\ requestedCut' = requestedCut
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = {}
    /\ compactionPlanning' = FALSE
    /\ compactionPending' = compactionPending \/ recoveryPending
    /\ cpuPressure' = cpuPressure
    /\ gpuPressure' = gpuPressure
    /\ recoveryPending' = recoveryPending
    /\ UNCHANGED <<
        quiescent, interactionState, workState, budgetState,
        staleAdmissionCount
       >>

RefreshSource ==
    /\ ~quiescent
    /\ sourceActive
    /\ generation < MaxEpoch
    /\ inventoryEpoch < MaxEpoch
    /\ routingId < MaxEpoch
    /\ sourceActive' = TRUE
    /\ viewEpoch' = viewEpoch
    /\ policyEpoch' = policyEpoch
    /\ inventoryEpoch' = inventoryEpoch + 1
    /\ routingId' = routingId + 1
    /\ generation' = generation + 1
    /\ covered' = {}
    /\ activeCut' = [o \in Objects |-> -1]
    /\ residentCut' = [o \in Objects |-> -1]
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = [o \in Objects |-> -1]
    /\ allocationView' = 0
    /\ allocationPolicy' = 0
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = {}
    /\ publicationPending' = TRUE
    /\ tasks' = tasks
    /\ results' = results
    /\ framePending' = framePending
    /\ compactionPlanning' = FALSE
    /\ compactionPending' = compactionPending \/ recoveryPending
    /\ cpuPressure' = cpuPressure
    /\ gpuPressure' = gpuPressure
    /\ recoveryPending' = recoveryPending
    /\ UNCHANGED <<
        quiescent, interactionState, budgetState, staleAdmissionCount
       >>

CancelGeneration ==
    /\ ~quiescent
    /\ sourceActive
    /\ sourceActive' = FALSE
    /\ viewEpoch' = viewEpoch
    /\ policyEpoch' = policyEpoch
    /\ inventoryEpoch' = inventoryEpoch
    /\ routingId' = routingId
    /\ generation' = generation
    /\ interactive' = FALSE
    /\ gestureActive' = FALSE
    /\ quietPending' = FALSE
    /\ covered' = {}
    /\ activeCut' = [o \in Objects |-> -1]
    /\ residentCut' = [o \in Objects |-> -1]
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = [o \in Objects |-> -1]
    /\ allocationView' = 0
    /\ allocationPolicy' = 0
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = {}
    /\ tasks' = tasks
    /\ results' = results
    /\ publicationPending' = TRUE
    /\ framePending' = framePending
    /\ renderBudget' = renderBudget
    /\ capacity' = capacity
    /\ calibrationPending' = FALSE
    /\ budgetRecoveryPending' = FALSE
    /\ cpuPressure' = FALSE
    /\ gpuPressure' = FALSE
    /\ recoveryPending' = FALSE
    /\ compactionPending' = FALSE
    /\ compactionPlanning' = FALSE
    /\ UNCHANGED <<quiescent, staleAdmissionCount>>

RestartSource ==
    /\ ~quiescent
    /\ ~sourceActive
    /\ generation < MaxEpoch
    /\ inventoryEpoch < MaxEpoch
    /\ routingId < MaxEpoch
    /\ sourceActive' = TRUE
    /\ viewEpoch' = viewEpoch
    /\ policyEpoch' = policyEpoch
    /\ inventoryEpoch' = inventoryEpoch + 1
    /\ routingId' = routingId + 1
    /\ generation' = generation + 1
    /\ covered' = {}
    /\ activeCut' = [o \in Objects |-> -1]
    /\ residentCut' = [o \in Objects |-> -1]
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = [o \in Objects |-> -1]
    /\ allocationView' = 0
    /\ allocationPolicy' = 0
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = {}
    /\ publicationPending' = TRUE
    /\ tasks' = tasks
    /\ results' = results
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, interactionState, budgetState, pressureState,
        staleAdmissionCount
       >>

BeginInteraction ==
    /\ ~quiescent
    /\ sourceActive
    /\ ~gestureActive
    /\ interactive' = TRUE
    /\ gestureActive' = TRUE
    /\ quietPending' = FALSE
    /\ UNCHANGED <<
        quiescent, identityState, geometryState, workState, budgetState,
        pressureState, staleAdmissionCount
       >>

EndGesture ==
    /\ gestureActive
    /\ interactive' = TRUE
    /\ gestureActive' = FALSE
    /\ quietPending' = TRUE
    /\ UNCHANGED <<
        quiescent, identityState, geometryState, workState, budgetState,
        pressureState, staleAdmissionCount
       >>

FinishQuiet ==
    /\ quietPending
    /\ interactive' = FALSE
    /\ gestureActive' = FALSE
    /\ quietPending' = FALSE
    /\ UNCHANGED <<
        quiescent, identityState, geometryState, workState, budgetState,
        pressureState, staleAdmissionCount
       >>

ChangeScanOrder ==
    /\ ~quiescent
    /\ \E order \in Orders:
        /\ scanOrder' = order
        /\ covered' = covered
        /\ activeCut' = activeCut
        /\ residentCut' = residentCut
        /\ requestedCut' = requestedCut
        /\ allocatedCut' = allocatedCut
        /\ allocationView' = allocationView
        /\ allocationPolicy' = allocationPolicy
        /\ memoryLimited' = memoryLimited
        /\ UNCHANGED <<
            quiescent, identityState, interactionState, workState,
            budgetState, pressureState, staleAdmissionCount
           >>

SubmitMinimum(o) ==
    /\ sourceActive
    /\ o \in Objects \ covered
    /\ ~tasks[o].valid
    /\ ~results[o].valid
    /\ tasks' = [tasks EXCEPT ![o] =
            MakeRequest(generation, viewEpoch, policyEpoch,
                inventoryEpoch, routingId, 0)]
    /\ results' = results
    /\ publicationPending' = publicationPending
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        budgetState, pressureState, staleAdmissionCount
       >>

SubmitAllocated(o) ==
    /\ sourceActive
    /\ covered = Objects
    /\ AllocationCurrent
    /\ o \notin memoryLimited
    /\ residentCut[o] < allocatedCut[o]
    /\ ~tasks[o].valid
    /\ ~results[o].valid
    /\ tasks' = [tasks EXCEPT ![o] =
            MakeRequest(generation, viewEpoch, policyEpoch,
                inventoryEpoch, routingId, allocatedCut[o])]
    /\ results' = results
    /\ publicationPending' = publicationPending
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        budgetState, pressureState, staleAdmissionCount
       >>

CompleteTask(o) ==
    /\ tasks[o].valid
    /\ ~results[o].valid
    /\ results' = [results EXCEPT ![o] = tasks[o]]
    /\ tasks' = [tasks EXCEPT ![o] = NoRequest]
    /\ publicationPending' = publicationPending
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        budgetState, pressureState, staleAdmissionCount
       >>

AdmitCurrentResult(o) ==
    /\ results[o].valid
    /\ RequestCurrent(results[o])
    /\ CanAdmit(o, results[o].cut)
    /\ residentCut' = ResidentWith(o, results[o].cut)
    /\ activeCut' = IF o \in covered THEN activeCut
            ELSE [activeCut EXCEPT ![o] = MaxNat(@, 0)]
    /\ covered' = covered \cup {o}
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = allocationView
    /\ allocationPolicy' = allocationPolicy
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = memoryLimited \ {o}
    /\ results' = [results EXCEPT ![o] = NoRequest]
    /\ tasks' = tasks
    /\ publicationPending' = publicationPending \/ o \notin covered
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, budgetState,
        pressureState, staleAdmissionCount
       >>

MarkMemoryLimited(o) ==
    /\ results[o].valid
    /\ RequestCurrent(results[o])
    /\ results[o].cut > 0
    /\ ~CanAdmit(o, results[o].cut)
    /\ memoryLimited' = memoryLimited \cup {o}
    /\ covered' = covered
    /\ activeCut' = activeCut
    /\ residentCut' = residentCut
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = allocationView
    /\ allocationPolicy' = allocationPolicy
    /\ scanOrder' = scanOrder
    /\ results' = [results EXCEPT ![o] = NoRequest]
    /\ tasks' = tasks
    /\ publicationPending' = publicationPending
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, budgetState,
        pressureState, staleAdmissionCount
       >>

RejectStaleResult(o) ==
    /\ results[o].valid
    /\ ~RequestCurrent(results[o])
    /\ results' = [results EXCEPT ![o] = NoRequest]
    /\ tasks' = tasks
    /\ covered' = IF InjectStaleAdmission THEN covered \cup {o} ELSE covered
    /\ residentCut' = IF InjectStaleAdmission
            THEN ResidentWith(o, results[o].cut) ELSE residentCut
    /\ activeCut' = IF InjectStaleAdmission
            THEN [activeCut EXCEPT ![o] = MaxNat(@, results[o].cut)]
            ELSE activeCut
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = allocationView
    /\ allocationPolicy' = allocationPolicy
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = memoryLimited
    /\ publicationPending' = publicationPending \/ InjectStaleAdmission
    /\ framePending' = framePending
    /\ staleAdmissionCount' = staleAdmissionCount +
            (IF InjectStaleAdmission THEN 1 ELSE 0)
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, budgetState, pressureState
       >>

AllocateScene ==
    /\ sourceActive
    /\ covered = Objects
    /\ ~interactive
    /\ ~AllocationCurrent
    /\ \E allocation \in BestAllocations:
        /\ allocatedCut' = allocation
        /\ allocationView' = viewEpoch
        /\ allocationPolicy' = policyEpoch
        /\ covered' = covered
        /\ activeCut' = activeCut
        /\ residentCut' = residentCut
        /\ requestedCut' = requestedCut
        /\ scanOrder' = scanOrder
        /\ memoryLimited' = {}
        /\ UNCHANGED <<
            quiescent, identityState, interactionState, workState,
            budgetState, pressureState, staleAdmissionCount
           >>

ApplyAllocation(o) ==
    /\ AllocationCurrent
    /\ o \notin memoryLimited
    /\ activeCut[o] # allocatedCut[o]
    /\ residentCut[o] >= allocatedCut[o]
    /\ activeCut' = [activeCut EXCEPT ![o] = allocatedCut[o]]
    /\ covered' = covered
    /\ residentCut' = residentCut
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = allocationView
    /\ allocationPolicy' = allocationPolicy
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = memoryLimited
    /\ publicationPending' = TRUE
    /\ tasks' = tasks
    /\ results' = results
    /\ framePending' = framePending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, budgetState,
        pressureState, staleAdmissionCount
       >>

RequestFrame ==
    /\ publicationPending
    /\ ~framePending
    /\ publicationPending' = FALSE
    /\ framePending' = TRUE
    /\ tasks' = tasks
    /\ results' = results
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        budgetState, pressureState, staleAdmissionCount
       >>

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ publicationPending' = publicationPending
    /\ tasks' = tasks
    /\ results' = results
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        budgetState, pressureState, staleAdmissionCount
       >>

ObservePressure ==
    /\ ~quiescent
    /\ sourceActive
    /\ \E cpu \in BOOLEAN, gpu \in BOOLEAN:
        /\ cpu \/ gpu
        /\ cpuPressure' = cpu
        /\ gpuPressure' = gpu
        /\ recoveryPending' = TRUE
        /\ compactionPending' = TRUE
        /\ compactionPlanning' = compactionPlanning
        /\ UNCHANGED <<
            quiescent, identityState, interactionState, geometryState,
            workState, budgetState, staleAdmissionCount
           >>

ClearPressure ==
    /\ ~quiescent
    /\ cpuPressure \/ gpuPressure
    /\ cpuPressure' = FALSE
    /\ gpuPressure' = FALSE
    /\ recoveryPending' = FALSE
    /\ compactionPending' = IF compactionPlanning THEN compactionPending
            ELSE FALSE
    /\ compactionPlanning' = compactionPlanning
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        workState, budgetState, staleAdmissionCount
       >>

ObserveCapacity ==
    /\ ~quiescent
    /\ sourceActive
    /\ \E observed \in CapacityStates:
        /\ capacity' = observed
        /\ calibrationPending' = TRUE
        /\ renderBudget' = renderBudget
        /\ budgetRecoveryPending' = budgetRecoveryPending
        /\ UNCHANGED <<
            quiescent, identityState, interactionState, geometryState,
            workState, pressureState, staleAdmissionCount
           >>

Calibrate ==
    /\ calibrationPending
    /\ LET nextBudget ==
            CASE capacity = "Overloaded" -> MaxNat(3, renderBudget - 1)
              [] capacity = "Headroom" ->
                    MinNat(MaxRenderBudget, renderBudget + 1)
              [] OTHER -> renderBudget
           changed == nextBudget # renderBudget
       IN
        /\ renderBudget' = nextBudget
        /\ capacity' = capacity
        /\ calibrationPending' = FALSE
        /\ budgetRecoveryPending' = budgetRecoveryPending \/
                (nextBudget < renderBudget
                 /\ TotalRender(activeCut) > EffectiveBudgetFor(nextBudget))
        /\ covered' = covered
        /\ activeCut' = activeCut
        /\ residentCut' = residentCut
        /\ requestedCut' = requestedCut
        /\ allocatedCut' = allocatedCut
        /\ allocationView' = IF changed THEN 0 ELSE allocationView
        /\ allocationPolicy' = IF changed THEN 0 ELSE allocationPolicy
        /\ scanOrder' = scanOrder
        /\ memoryLimited' = IF changed THEN {} ELSE memoryLimited
        /\ UNCHANGED <<
            quiescent, identityState, interactionState, workState,
            pressureState, staleAdmissionCount
           >>

FinishBudgetRecovery ==
    /\ budgetRecoveryPending
    /\ AllocationCurrent
    /\ CurrentVisualSettled
    /\ TotalRender(activeCut) <= EffectiveRenderBudget
    /\ budgetRecoveryPending' = FALSE
    /\ renderBudget' = renderBudget
    /\ capacity' = capacity
    /\ calibrationPending' = calibrationPending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        workState, pressureState, staleAdmissionCount
       >>

CompactionPrerequisites ==
    /\ compactionPending
    /\ sourceActive
    /\ covered = Objects
    /\ ~interactive
    /\ ~gestureActive
    /\ ~quietPending
    /\ ~HasTask
    /\ ~HasResult
    /\ ~publicationPending
    /\ ~framePending
    /\ ~calibrationPending
    /\ ~budgetRecoveryPending
    /\ AllocationCurrent
    /\ CurrentVisualSettled

PlanCompaction ==
    /\ CompactionPrerequisites
    /\ ~compactionPlanning
    /\ compactionPlanning' = TRUE
    /\ cpuPressure' = cpuPressure
    /\ gpuPressure' = gpuPressure
    /\ recoveryPending' = recoveryPending
    /\ compactionPending' = compactionPending
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, geometryState,
        workState, budgetState, staleAdmissionCount
       >>

FinishCompaction ==
    /\ compactionPlanning
    /\ residentCut' = [o \in Objects |-> activeCut[o]]
    /\ covered' = covered
    /\ activeCut' = activeCut
    /\ requestedCut' = requestedCut
    /\ allocatedCut' = allocatedCut
    /\ allocationView' = allocationView
    /\ allocationPolicy' = allocationPolicy
    /\ scanOrder' = scanOrder
    /\ memoryLimited' = memoryLimited
    /\ cpuPressure' = cpuPressure
    /\ gpuPressure' = gpuPressure
    /\ recoveryPending' = FALSE
    /\ compactionPending' = FALSE
    /\ compactionPlanning' = FALSE
    /\ UNCHANGED <<
        quiescent, identityState, interactionState, workState, budgetState,
        staleAdmissionCount
       >>

EnvironmentAction ==
    InvalidateView \/ ChangePolicy \/ RefreshSource \/ CancelGeneration
    \/ RestartSource \/ BeginInteraction \/ ChangeScanOrder
    \/ ObservePressure \/ ClearPressure \/ ObserveCapacity

InternalAction ==
    EndGesture \/ FinishQuiet
    \/ (\E o \in Objects:
        SubmitMinimum(o) \/ SubmitAllocated(o) \/ CompleteTask(o)
        \/ AdmitCurrentResult(o) \/ MarkMemoryLimited(o)
        \/ RejectStaleResult(o) \/ ApplyAllocation(o))
    \/ AllocateScene \/ RequestFrame \/ CompleteFrame \/ Calibrate
    \/ FinishBudgetRecovery \/ PlanCompaction \/ FinishCompaction

Idle ==
    /\ quiescent
    /\ Terminal
    /\ UNCHANGED vars

Next == QuiesceEnvironment \/ EnvironmentAction \/ InternalAction \/ Idle

Spec == Init /\ [][Next]_vars

LiveSpec ==
    Spec
    /\ WF_vars(QuiesceEnvironment)
    /\ WF_vars(EndGesture)
    /\ WF_vars(FinishQuiet)
    /\ WF_vars(AllocateScene)
    /\ WF_vars(RequestFrame)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(Calibrate)
    /\ WF_vars(FinishBudgetRecovery)
    /\ WF_vars(PlanCompaction)
    /\ WF_vars(FinishCompaction)
    /\ \A o \in Objects:
        /\ WF_vars(SubmitMinimum(o))
        /\ WF_vars(SubmitAllocated(o))
        /\ WF_vars(CompleteTask(o))
        /\ WF_vars(AdmitCurrentResult(o))
        /\ WF_vars(MarkMemoryLimited(o))
        /\ WF_vars(RejectStaleResult(o))
        /\ WF_vars(ApplyAllocation(o))

TypeOK ==
    /\ quiescent \in BOOLEAN
    /\ sourceActive \in BOOLEAN
    /\ viewEpoch \in 1..MaxEpoch
    /\ policyEpoch \in 1..MaxEpoch
    /\ inventoryEpoch \in 1..MaxEpoch
    /\ routingId \in 1..MaxEpoch
    /\ generation \in 1..MaxEpoch
    /\ interactive \in BOOLEAN
    /\ gestureActive \in BOOLEAN
    /\ quietPending \in BOOLEAN
    /\ covered \subseteq Objects
    /\ activeCut \in [Objects -> CutDomain]
    /\ residentCut \in [Objects -> CutDomain]
    /\ requestedCut \in [Objects -> Cuts]
    /\ allocatedCut \in [Objects -> CutDomain]
    /\ allocationView \in 0..MaxEpoch
    /\ allocationPolicy \in 0..MaxEpoch
    /\ scanOrder \in Orders
    /\ memoryLimited \subseteq Objects
    /\ \A o \in Objects: RequestType(tasks[o]) /\ RequestType(results[o])
    /\ publicationPending \in BOOLEAN
    /\ framePending \in BOOLEAN
    /\ renderBudget \in 3..MaxRenderBudget
    /\ capacity \in CapacityStates
    /\ calibrationPending \in BOOLEAN
    /\ budgetRecoveryPending \in BOOLEAN
    /\ cpuPressure \in BOOLEAN
    /\ gpuPressure \in BOOLEAN
    /\ recoveryPending \in BOOLEAN
    /\ compactionPending \in BOOLEAN
    /\ compactionPlanning \in BOOLEAN
    /\ staleAdmissionCount \in Nat

ActiveResidentOrder ==
    \A o \in Objects: activeCut[o] <= residentCut[o]

CoverageIsDrawable ==
    \A o \in covered: activeCut[o] >= 0 /\ residentCut[o] >= activeCut[o]

NoStaleAdmissions == staleAdmissionCount = 0

IrreducibleResidentOverage ==
    TotalMemory(residentCut) > ResidentLimit
    /\ \A o \in Objects: residentCut[o] <= 0

ResidentBudgetSafety ==
    TotalMemory(residentCut) <= ResidentLimit \/ IrreducibleResidentOverage

ActiveBudgetSafety ==
    TotalRender(activeCut) <= EffectiveRenderBudget \/ budgetRecoveryPending

AllocationEpochSafety ==
    AllocationCurrent => allocatedCut \in BestAllocations

ProtectedFloorRespected ==
    AllocationCurrent /\ ProtectedFloorFeasible =>
        \A o \in Objects: allocatedCut[o] >= ProtectedCut[o]

ImportanceOrdering ==
    AllocationCurrent =>
        allocatedCut["A"] >= allocatedCut["B"]
        /\ allocatedCut["B"] >= allocatedCut["C"]

LocallyMaximalAllocation ==
    AllocationCurrent =>
        \A o \in Objects:
            allocatedCut[o] < requestedCut[o] =>
                TotalRender(allocatedCut)
                - RenderAt(allocatedCut[o])
                + RenderAt(allocatedCut[o] + 1) > EffectiveRenderBudget

CompactionSafety ==
    compactionPlanning =>
        sourceActive /\ covered = Objects
        /\ \A o \in Objects: activeCut[o] <= residentCut[o]

StableHasNoForeground ==
    ReportedStable => ForegroundWitnesses = {}

StableOutcomeClassified ==
    ReportedStable =>
        FullRequestedQuality \/ PerformanceLimited \/ MemoryLimitedTerminal

EveryNonterminalHasWitness ==
    ~Terminal => ForegroundWitnesses # {} \/ compactionPending

QuiescentTermination == quiescent ~> Terminal

=============================================================================
