----------------------- MODULE ObolProgressivePipeline --------------------
\* Canonical control-plane contract for progressive CAD presentation.
\*
\* Workload names are model-checking profiles, not production modes.  The
\* same discovery, availability, planning, presentation, and convergence
\* protocol handles all of them.  Two occurrences/assets are a symmetry
\* reduction for cardinalities which range from one to hundreds of thousands.
\*
\* This model deliberately excludes triangle construction, numeric visual
\* error, and elapsed time.  Those belong to pure planner tests, renderer
\* boundary tests, perf, and image qualification.  It does prove that those
\* mechanisms cannot create a second control owner, conceal an error as
\* readiness, starve one of several large assets, or balance the same evidence
\* repeatedly.

EXTENDS Naturals, FiniteSets, TLC

CONSTANT MaxViewEpoch, MaxPolicyEpoch

Assets == {"assetA", "assetB"}
Occurrences == {"hero", "peer"}
Workloads == {"fewSmall", "manySmall", "singleLarge",
              "multiLarge", "sharedLarge"}
AvailabilityStates == {"unknown", "minimum", "rich", "failed"}
Representations == {"none", "box", "aggregate", "preview", "mesh"}
MeshKinds == {"none", "direct", "pop"}
WorkGoals == {"none", "minimum", "rich"}
TerminalOutcomes == {"ready", "constrained", "error"}
InitialBudgets == 0..3

RevisionStamp(inventory, available, view, policy, capacity) ==
    [inventory |-> inventory,
     availability |-> available,
     view |-> view,
     policy |-> policy,
     capacity |-> capacity]

ActiveOccurrences(workload) ==
    IF workload = "singleLarge" THEN {"hero"} ELSE Occurrences

ActiveAssets(workload) ==
    IF workload \in {"singleLarge", "sharedLarge"}
    THEN {"assetA"}
    ELSE Assets

AssetFor(workload, occurrence) ==
    IF workload = "sharedLarge" THEN "assetA"
    ELSE IF occurrence = "hero" THEN "assetA" ELSE "assetB"

LargeAsset(workload, asset) ==
    workload \in {"singleLarge", "multiLarge", "sharedLarge"} /\
    asset \in ActiveAssets(workload)

DirectSafe(workload) == workload = "fewSmall"

Subpixel(workload, occurrence) ==
    workload = "manySmall" /\ occurrence = "peer"

Cost(workload, occurrence) ==
    IF LargeAsset(workload, AssetFor(workload, occurrence))
    THEN 2 ELSE 1

RepRank(rep) ==
    IF rep = "none" THEN 0
    ELSE IF rep = "box" THEN 1
    ELSE IF rep \in {"aggregate", "preview"} THEN 2
    ELSE 3

VARIABLES workload,
          viewEpoch,
          inputOpen,
          discovered,
          profileCertified,
          availability,
          enrichmentExhausted,
          workingAsset,
          workGoal,
          committedRep,
          committedQuality,
          committedKind,
          candidateRep,
          candidateQuality,
          candidateKind,
          framePending,
          overviewVisible,
          budget,
          capacityMeasured,
          inventoryRevision,
          availabilityRevision,
          policyRevision,
          capacityRevision,
          plannedRevisionStamp,
          planningRevisionStamp,
          planning,
          considered,
          remainingBudget,
          lastCommitEpoch,
          previousCommitEpoch,
          previousCommittedRep,
          previousCommittedQuality

vars == <<workload, viewEpoch, inputOpen, discovered, profileCertified,
          availability, enrichmentExhausted, workingAsset, workGoal,
          committedRep, committedQuality, committedKind,
          candidateRep, candidateQuality, candidateKind, framePending,
          overviewVisible, budget, capacityMeasured, inventoryRevision,
          availabilityRevision, policyRevision, capacityRevision,
          plannedRevisionStamp, planningRevisionStamp, planning, considered,
          remainingBudget,
          lastCommitEpoch, previousCommitEpoch, previousCommittedRep,
          previousCommittedQuality>>

AllDiscovered ==
    \A occurrence \in ActiveOccurrences(workload): discovered[occurrence]

AvailabilityTerminal(asset) ==
    availability[asset] \in {"rich", "failed"} \/
    (availability[asset] = "minimum" /\ enrichmentExhausted[asset])

AllAvailabilityTerminal ==
    \A asset \in ActiveAssets(workload): AvailabilityTerminal(asset)

HasTerminalFailure ==
    \E occurrence \in ActiveOccurrences(workload):
        availability[AssetFor(workload, occurrence)] = "failed" /\
        committedRep[occurrence] \in {"none", "box"}

HasStructuralBox ==
    \E occurrence \in ActiveOccurrences(workload):
        committedRep[occurrence] \in {"none", "box"}

HasQualityDebt ==
    \E occurrence \in ActiveOccurrences(workload):
        /\ ~Subpixel(workload, occurrence)
        /\ (committedRep[occurrence] # "mesh" \/
             committedQuality[occurrence] < 2)

CurrentRevisionStamp ==
    RevisionStamp(inventoryRevision, availabilityRevision, viewEpoch,
                  policyRevision, capacityRevision)

PipelineBusy ==
    inputOpen \/ ~AllDiscovered \/ ~profileCertified \/
    ~AllAvailabilityTerminal \/ workingAsset # "none" \/
    planning \/ framePending \/
    plannedRevisionStamp # CurrentRevisionStamp

Outcome ==
    IF PipelineBusy THEN "active"
    ELSE IF HasTerminalFailure THEN "error"
    ELSE IF HasStructuralBox THEN "error"
    ELSE IF HasQualityDebt THEN "constrained"
    ELSE "ready"

TypeOK ==
    /\ workload \in Workloads
    /\ viewEpoch \in 0..MaxViewEpoch
    /\ inputOpen \in BOOLEAN
    /\ discovered \in [Occurrences -> BOOLEAN]
    /\ profileCertified \in BOOLEAN
    /\ availability \in [Assets -> AvailabilityStates]
    /\ enrichmentExhausted \in [Assets -> BOOLEAN]
    /\ workingAsset \in Assets \union {"none"}
    /\ workGoal \in WorkGoals
    /\ committedRep \in [Occurrences -> Representations]
    /\ committedQuality \in [Occurrences -> 0..2]
    /\ committedKind \in [Occurrences -> MeshKinds]
    /\ candidateRep \in [Occurrences -> Representations]
    /\ candidateQuality \in [Occurrences -> 0..2]
    /\ candidateKind \in [Occurrences -> MeshKinds]
    /\ framePending \in BOOLEAN
    /\ overviewVisible \in BOOLEAN
    /\ budget \in 0..3
    /\ capacityMeasured \in BOOLEAN
    /\ inventoryRevision \in Nat
    /\ availabilityRevision \in Nat
    /\ policyRevision \in 0..MaxPolicyEpoch
    /\ capacityRevision \in Nat
    /\ plannedRevisionStamp \in
         [inventory: Nat, availability: Nat, view: Nat,
          policy: Nat, capacity: Nat]
    /\ planningRevisionStamp \in
         [inventory: Nat, availability: Nat, view: Nat,
          policy: Nat, capacity: Nat]
    /\ planning \in BOOLEAN
    /\ considered \subseteq Occurrences
    /\ remainingBudget \in 0..3
    /\ lastCommitEpoch \in 0..MaxViewEpoch
    /\ previousCommitEpoch \in 0..MaxViewEpoch
    /\ previousCommittedRep \in [Occurrences -> Representations]
    /\ previousCommittedQuality \in [Occurrences -> 0..2]
    /\ Outcome \in TerminalOutcomes \union {"active"}

Init ==
    /\ workload \in Workloads
    /\ viewEpoch = 0
    /\ inputOpen = FALSE
    /\ discovered = [occurrence \in Occurrences |-> FALSE]
    /\ profileCertified = FALSE
    /\ availability = [asset \in Assets |->
           IF asset \in ActiveAssets(workload) THEN "unknown" ELSE "rich"]
    /\ enrichmentExhausted = [asset \in Assets |-> FALSE]
    /\ workingAsset = "none"
    /\ workGoal = "none"
    /\ committedRep = [occurrence \in Occurrences |-> "none"]
    /\ committedQuality = [occurrence \in Occurrences |-> 0]
    /\ committedKind = [occurrence \in Occurrences |-> "none"]
    /\ candidateRep = committedRep
    /\ candidateQuality = committedQuality
    /\ candidateKind = committedKind
    /\ framePending = FALSE
    /\ overviewVisible = TRUE
    /\ budget \in InitialBudgets
    /\ capacityMeasured = FALSE
    /\ inventoryRevision = 1
    /\ availabilityRevision = 1
    /\ policyRevision = 0
    /\ capacityRevision = 1
    /\ plannedRevisionStamp = RevisionStamp(0, 0, 0, 0, 0)
    /\ planningRevisionStamp = RevisionStamp(0, 0, 0, 0, 0)
    /\ planning = FALSE
    /\ considered = {}
    /\ remainingBudget = budget
    /\ lastCommitEpoch = 0
    /\ previousCommitEpoch = 0
    /\ previousCommittedRep = committedRep
    /\ previousCommittedQuality = committedQuality

Discover(occurrence) ==
    /\ occurrence \in ActiveOccurrences(workload)
    /\ ~discovered[occurrence]
    /\ ~framePending
    /\ ~planning
    /\ discovered' = [discovered EXCEPT ![occurrence] = TRUE]
    /\ candidateRep' = [committedRep EXCEPT ![occurrence] = "box"]
    /\ candidateQuality' = committedQuality
    /\ candidateKind' = committedKind
    /\ framePending' = TRUE
    /\ inventoryRevision' = inventoryRevision + 1
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, profileCertified,
                    availability, enrichmentExhausted, workingAsset,
                    workGoal, committedRep, committedQuality, committedKind,
                    overviewVisible, budget, capacityMeasured,
                    availabilityRevision, policyRevision, capacityRevision,
                    plannedRevisionStamp, planningRevisionStamp, planning,
                    considered, remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

CertifyProfile ==
    /\ AllDiscovered
    /\ ~profileCertified
    /\ ~framePending
    /\ profileCertified' = TRUE
    /\ inventoryRevision' = inventoryRevision + 1
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    availability, enrichmentExhausted, workingAsset,
                    workGoal, committedRep, committedQuality, committedKind,
                    candidateRep, candidateQuality, candidateKind,
                    framePending, overviewVisible, budget, capacityMeasured,
                    availabilityRevision, policyRevision, capacityRevision,
                    plannedRevisionStamp, planningRevisionStamp, planning,
                    considered,
                    remainingBudget, lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

StartMinimum(asset) ==
    /\ profileCertified
    /\ asset \in ActiveAssets(workload)
    /\ availability[asset] = "unknown"
    /\ workingAsset = "none"
    /\ workingAsset' = asset
    /\ workGoal' = "minimum"
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    committedRep, committedQuality, committedKind,
                    candidateRep, candidateQuality, candidateKind,
                    framePending, overviewVisible, budget, capacityMeasured,
                    inventoryRevision, availabilityRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

FinishMinimum ==
    /\ workingAsset \in Assets
    /\ workGoal = "minimum"
    /\ \/ /\ availability' = [availability EXCEPT
                    ![workingAsset] = IF DirectSafe(workload)
                                      THEN "rich" ELSE "minimum"]
           /\ UNCHANGED enrichmentExhausted
       \/ /\ availability' = [availability EXCEPT
                    ![workingAsset] = "failed"]
           /\ enrichmentExhausted' = [enrichmentExhausted EXCEPT
                    ![workingAsset] = TRUE]
    /\ workingAsset' = "none"
    /\ workGoal' = "none"
    /\ availabilityRevision' = availabilityRevision + 1
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, committedRep, committedQuality,
                    committedKind, candidateRep, candidateQuality,
                    candidateKind, framePending, overviewVisible, budget,
                    capacityMeasured, inventoryRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

StartRich(asset) ==
    /\ profileCertified
    /\ asset \in ActiveAssets(workload)
    /\ availability[asset] = "minimum"
    /\ ~enrichmentExhausted[asset]
    /\ workingAsset = "none"
    /\ workingAsset' = asset
    /\ workGoal' = "rich"
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    committedRep, committedQuality, committedKind,
                    candidateRep, candidateQuality, candidateKind,
                    framePending, overviewVisible, budget, capacityMeasured,
                    inventoryRevision, availabilityRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

FinishRich ==
    /\ workingAsset \in Assets
    /\ workGoal = "rich"
    /\ \/ /\ availability' = [availability EXCEPT
                    ![workingAsset] = "rich"]
           /\ UNCHANGED enrichmentExhausted
       \/ /\ UNCHANGED availability
           /\ enrichmentExhausted' = [enrichmentExhausted EXCEPT
                    ![workingAsset] = TRUE]
    /\ workingAsset' = "none"
    /\ workGoal' = "none"
    /\ availabilityRevision' = availabilityRevision + 1
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, committedRep, committedQuality,
                    committedKind, candidateRep, candidateQuality,
                    candidateKind, framePending, overviewVisible, budget,
                    capacityMeasured, inventoryRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

StartPlan ==
    /\ ~inputOpen
    /\ profileCertified
    /\ AllDiscovered
    /\ ~planning
    /\ ~framePending
    /\ plannedRevisionStamp # CurrentRevisionStamp
    /\ planning' = TRUE
    /\ planningRevisionStamp' = CurrentRevisionStamp
    /\ considered' = Occurrences \ ActiveOccurrences(workload)
    /\ remainingBudget' = budget
    /\ candidateRep' = committedRep
    /\ candidateQuality' = committedQuality
    /\ candidateKind' = committedKind
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, framePending, overviewVisible, budget,
                    capacityMeasured, inventoryRevision,
                    availabilityRevision, policyRevision, capacityRevision,
                    plannedRevisionStamp, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

\* A camera edge changes demand; it does not erase a complete retained
\* presentation.  The next stable plan starts from that presentation and may
\* improve it.  A future model which admits quality regression must first add
\* a typed memory/deadline constraint revision and make that proof explicit.
CanKeepRetainedMesh(occurrence) ==
    committedRep[occurrence] = "mesh"

ChosenRep(occurrence) ==
    IF Subpixel(workload, occurrence) THEN "aggregate"
    ELSE IF availability[AssetFor(workload, occurrence)] = "failed"
         THEN committedRep[occurrence]
    ELSE IF availability[AssetFor(workload, occurrence)] \in
                {"minimum", "rich"}
         THEN IF CanKeepRetainedMesh(occurrence) \/
                     remainingBudget >= Cost(workload, occurrence)
              THEN "mesh" ELSE "preview"
    ELSE committedRep[occurrence]

ChosenQuality(occurrence) ==
    IF ChosenRep(occurrence) # "mesh" THEN 0
    ELSE IF DirectSafe(workload) \/
            availability[AssetFor(workload, occurrence)] = "rich"
         THEN 2
    ELSE IF CanKeepRetainedMesh(occurrence)
         THEN committedQuality[occurrence]
    ELSE 1

ChosenKind(occurrence) ==
    IF ChosenRep(occurrence) # "mesh" THEN "none"
    ELSE IF DirectSafe(workload) THEN "direct" ELSE "pop"

PlanOccurrence(occurrence) ==
    /\ planning
    /\ planningRevisionStamp = CurrentRevisionStamp
    /\ occurrence \in ActiveOccurrences(workload)
    /\ occurrence \notin considered
    /\ occurrence = "hero" \/ "hero" \in considered
    /\ candidateRep' = [candidateRep EXCEPT
            ![occurrence] = ChosenRep(occurrence)]
    /\ candidateQuality' = [candidateQuality EXCEPT
            ![occurrence] = ChosenQuality(occurrence)]
    /\ candidateKind' = [candidateKind EXCEPT
            ![occurrence] = ChosenKind(occurrence)]
    /\ remainingBudget' =
         IF ChosenRep(occurrence) = "mesh" /\
                ~CanKeepRetainedMesh(occurrence)
         THEN remainingBudget - Cost(workload, occurrence)
         ELSE remainingBudget
    /\ considered' = considered \union {occurrence}
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, framePending, overviewVisible, budget,
                    capacityMeasured, inventoryRevision,
                    availabilityRevision, policyRevision, capacityRevision,
                    plannedRevisionStamp, planningRevisionStamp, planning,
                    lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

FinishPlan ==
    /\ planning
    /\ planningRevisionStamp = CurrentRevisionStamp
    /\ considered = Occurrences
    /\ planning' = FALSE
    /\ plannedRevisionStamp' = planningRevisionStamp
    /\ framePending' =
         (candidateRep # committedRep \/
          candidateQuality # committedQuality \/
          candidateKind # committedKind)
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, candidateRep, candidateQuality,
                    candidateKind, overviewVisible, budget,
                    capacityMeasured, inventoryRevision,
                    availabilityRevision, policyRevision, capacityRevision,
                    planningRevisionStamp, considered,
                    remainingBudget, lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

AbortStalePlan ==
    /\ planning
    /\ planningRevisionStamp # CurrentRevisionStamp
    /\ planning' = FALSE
    /\ considered' = {}
    /\ remainingBudget' = budget
    /\ candidateRep' = committedRep
    /\ candidateQuality' = committedQuality
    /\ candidateKind' = committedKind
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, framePending, overviewVisible, budget,
                    capacityMeasured, inventoryRevision,
                    availabilityRevision, policyRevision, capacityRevision,
                    plannedRevisionStamp, planningRevisionStamp,
                    lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

CompleteFrame ==
    /\ framePending
    /\ committedRep' = candidateRep
    /\ committedQuality' = candidateQuality
    /\ committedKind' = candidateKind
    /\ framePending' = FALSE
    /\ overviewVisible' =
         ~\A occurrence \in ActiveOccurrences(workload):
             candidateRep[occurrence] # "none"
    /\ previousCommitEpoch' = lastCommitEpoch
    /\ previousCommittedRep' = committedRep
    /\ previousCommittedQuality' = committedQuality
    /\ lastCommitEpoch' = viewEpoch
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, candidateRep, candidateQuality,
                    candidateKind, budget, capacityMeasured,
                    inventoryRevision, availabilityRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget>>

MeasureCapacity(newBudget) ==
    /\ newBudget \in budget..3
    /\ ~inputOpen
    /\ profileCertified
    /\ AllAvailabilityTerminal
    /\ ~capacityMeasured
    /\ ~planning
    /\ ~framePending
    /\ budget' = newBudget
    /\ capacityMeasured' = TRUE
    /\ capacityRevision' = capacityRevision + 1
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, candidateRep, candidateQuality,
                    candidateKind, framePending, overviewVisible,
                    inventoryRevision, availabilityRevision, policyRevision,
                    plannedRevisionStamp, planningRevisionStamp, planning,
                    considered,
                    remainingBudget, lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

BeginInput ==
    /\ ~inputOpen
    /\ viewEpoch < MaxViewEpoch
    /\ ~framePending
    /\ viewEpoch' = viewEpoch + 1
    /\ inputOpen' = TRUE
    /\ capacityMeasured' = FALSE
    /\ UNCHANGED <<workload, discovered, profileCertified, availability,
                    enrichmentExhausted, workingAsset, workGoal,
                    committedRep, committedQuality, committedKind,
                    candidateRep, candidateQuality, candidateKind,
                    framePending, overviewVisible, budget,
                    inventoryRevision, availabilityRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

EndInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ UNCHANGED <<workload, viewEpoch, discovered, profileCertified,
                    availability, enrichmentExhausted, workingAsset,
                    workGoal, committedRep, committedQuality, committedKind,
                    candidateRep, candidateQuality, candidateKind,
                    framePending, overviewVisible, budget, capacityMeasured,
                    inventoryRevision, availabilityRevision, policyRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch,
                    previousCommitEpoch, previousCommittedRep,
                    previousCommittedQuality>>

ChangePolicy ==
    /\ ~inputOpen
    /\ ~framePending
    /\ policyRevision < MaxPolicyEpoch
    /\ policyRevision' = policyRevision + 1
    /\ capacityMeasured' = FALSE
    /\ UNCHANGED <<workload, viewEpoch, inputOpen, discovered,
                    profileCertified, availability, enrichmentExhausted,
                    workingAsset, workGoal, committedRep, committedQuality,
                    committedKind, candidateRep, candidateQuality,
                    candidateKind, framePending, overviewVisible, budget,
                    inventoryRevision, availabilityRevision,
                    capacityRevision, plannedRevisionStamp,
                    planningRevisionStamp, planning, considered,
                    remainingBudget, lastCommitEpoch, previousCommitEpoch,
                    previousCommittedRep, previousCommittedQuality>>

Next ==
    \/ \E occurrence \in Occurrences: Discover(occurrence)
    \/ CertifyProfile
    \/ \E asset \in Assets: StartMinimum(asset)
    \/ FinishMinimum
    \/ \E asset \in Assets: StartRich(asset)
    \/ FinishRich
    \/ StartPlan
    \/ \E occurrence \in Occurrences: PlanOccurrence(occurrence)
    \/ FinishPlan
    \/ AbortStalePlan
    \/ CompleteFrame
    \/ \E newBudget \in InitialBudgets: MeasureCapacity(newBudget)
    \/ BeginInput
    \/ EndInput
    \/ ChangePolicy

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ \A occurrence \in Occurrences: WF_vars(Discover(occurrence))
    /\ WF_vars(CertifyProfile)
    /\ \A asset \in Assets: WF_vars(StartMinimum(asset))
    /\ WF_vars(FinishMinimum)
    /\ \A asset \in Assets: WF_vars(StartRich(asset))
    /\ WF_vars(FinishRich)
    /\ WF_vars(StartPlan)
    /\ \A occurrence \in Occurrences: WF_vars(PlanOccurrence(occurrence))
    /\ WF_vars(FinishPlan)
    /\ WF_vars(AbortStalePlan)
    /\ WF_vars(CompleteFrame)
    /\ WF_vars(EndInput)

CommittedCoverageNeverDisappears ==
    overviewVisible \/
    \A occurrence \in ActiveOccurrences(workload):
        committedRep[occurrence] # "none"

PreparationMemoryBounded ==
    Cardinality({asset \in Assets: workingAsset = asset}) <= 1

SharedAssetPreparedOnce ==
    workload = "sharedLarge" => "assetB" \notin ActiveAssets(workload)

LargeWorkloadsHaveNoDirectMeshes ==
    workload # "fewSmall" =>
        \A occurrence \in ActiveOccurrences(workload):
            committedKind[occurrence] # "direct"

MeshKindIsExact ==
    \A occurrence \in Occurrences:
        (committedRep[occurrence] = "mesh") =
            (committedKind[occurrence] \in {"direct", "pop"})

TerminalSuccessHasNoStructuralBoxes ==
    Outcome \in {"ready", "constrained"} => ~HasStructuralBox

ErrorIsNeverReady == Outcome = "error" => Outcome # "ready"

PlanningRequiresNewEvidence ==
    planning => plannedRevisionStamp # planningRevisionStamp

NoRepeatedPlanWithoutEvidence ==
    ~planning /\ plannedRevisionStamp = CurrentRevisionStamp /\
    ~framePending => ~ENABLED StartPlan

StalePlanCannotAdvance ==
    planning /\ planningRevisionStamp # CurrentRevisionStamp =>
        /\ ~\E occurrence \in Occurrences: ENABLED PlanOccurrence(occurrence)
        /\ ~ENABLED FinishPlan
        /\ ENABLED AbortStalePlan

SameEpochFrameDoesNotRegress ==
    lastCommitEpoch = previousCommitEpoch =>
        \A occurrence \in ActiveOccurrences(workload):
            /\ RepRank(committedRep[occurrence]) >=
               RepRank(previousCommittedRep[occurrence])
            /\ committedQuality[occurrence] >=
               previousCommittedQuality[occurrence]

\* View and input epochs are demand revisions, not capacity proofs.  Existing
\* complete mesh data remains the starting point after pose, translation, or
\* zoom; ordinary convergence cannot replace it with a coarser stable commit.
\* A later extension may weaken this only behind an explicit typed capacity or
\* memory constraint carried by the committed transaction.
RetainedMeshDoesNotRegressAcrossViewEpochs ==
    \A occurrence \in ActiveOccurrences(workload):
        previousCommittedRep[occurrence] = "mesh" =>
            /\ committedRep[occurrence] = "mesh"
            /\ committedQuality[occurrence] >=
               previousCommittedQuality[occurrence]

ActivePipelineHasProgressWitness ==
    Outcome = "active" =>
        \/ ENABLED CertifyProfile
        \/ ENABLED FinishMinimum
        \/ ENABLED FinishRich
        \/ ENABLED StartPlan
        \/ ENABLED FinishPlan
        \/ ENABLED AbortStalePlan
        \/ ENABLED CompleteFrame
        \/ ENABLED EndInput
        \/ \E occurrence \in Occurrences:
              ENABLED Discover(occurrence) \/
              ENABLED PlanOccurrence(occurrence)
        \/ \E asset \in Assets:
              ENABLED StartMinimum(asset) \/ ENABLED StartRich(asset)
        \/ \E newBudget \in InitialBudgets:
              ENABLED MeasureCapacity(newBudget)

EventuallyTerminalAfterFinalInput ==
    []((viewEpoch = MaxViewEpoch /\ ~inputOpen) =>
       <>(Outcome \in TerminalOutcomes))

MultiLargeAssetsDoNotStarve ==
    []((workload = "multiLarge" /\ viewEpoch = MaxViewEpoch /\ ~inputOpen) =>
       <>(AvailabilityTerminal("assetA") /\
          AvailabilityTerminal("assetB")))

=============================================================================
