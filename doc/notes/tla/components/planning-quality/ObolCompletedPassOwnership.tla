----------------------- MODULE ObolCompletedPassOwnership ----------------------
\* Multi-pass successor transaction for a completed progressive-planning pass.
\*
\* Selection is made from one immutable completion snapshot before any owner
\* mutates evidence.  Resident growth and coverage have precedence over
\* structural-frontier repair, presentation handoff, and generic capacity
\* calibration.  Structural-only occurrences cannot participate in an
\* occurrence allocation or renderer-capacity sample; their exact projected
\* frontier therefore invalidates any stale sample and owns a finite repair
\* successor before either operation.  Every effect which
\* starts a successor pass clears the complete pass-annotation transaction.
\* A scene-wide allocation successor also consumes any selective source-delta
\* plan owned by the pass which just completed.
\* Capacity samples of the same discrete population update evidence without
\* advancing the semantic capacity revision.
\*
\* Production refinement:
\*  - BObolLodAvailabilityScheduler selects resident growth and coverage;
\*  - BObolLodPresentationPolicy::completedPassSelection selects handoff or
\*    capacity from a normalized immutable pass snapshot;
\*  - BObolLodRetainedPassAnnotations owns the pass-scoped annotations; and
\*  - BObolLodRevisionContract owns semantic revision advancement.

EXTENDS Naturals, TLC

Phases == {"completed", "execute", "terminal"}
Handoffs == {"none", "presentation", "allocation"}
Owners == {"none", "growth", "coverage", "structural", "presentation",
           "allocation", "capacity", "observation"}
Effects == {"none", "drainGrowth", "retryCoverage", "repairStructural",
            "presentHandoff",
            "allocateHandoff", "finishHandoff", "presentCut",
            "changePopulation", "recordSample", "retireObservation"}
Populations == 0..1
\* At most growth, coverage, three handoff effects (presentation, allocation,
\* finish), one structural repair, three capacity effects, and one retained
\* observation retirement can precede terminal selection in this
\* symmetry-reduced transaction.
MaximumCompletedPassCount == 10

VARIABLES phase,
          growthPending,
          coveragePending,
          structuralPending,
          retainedRefinementPending,
          selectiveDelta,
          changedCut,
          rescanAfterFrame,
          handoff,
          allocationCertified,
          ceilingFreeSample,
          capacityPending,
          capacityEligible,
          capacityProducer,
          capacityPopulation,
          capacityCandidate,
          initialCapacityPopulation,
          capacityRevision,
          owner,
          effect,
          completedPassCount

vars == <<phase, growthPending, coveragePending, structuralPending,
          retainedRefinementPending,
          selectiveDelta, changedCut,
          rescanAfterFrame, handoff, allocationCertified, ceilingFreeSample,
          capacityPending, capacityEligible, capacityProducer,
          capacityPopulation, capacityCandidate, initialCapacityPopulation,
          capacityRevision, owner, effect, completedPassCount>>

CleanPass == ~changedCut /\ ~rescanAfterFrame
NormalizedCleanPass == CleanPass \/ ceilingFreeSample

ExpectedOwner ==
    IF growthPending THEN "growth"
    ELSE IF coveragePending THEN "coverage"
    ELSE IF structuralPending THEN "structural"
    ELSE IF NormalizedCleanPass /\ handoff = "presentation" /\
            (~capacityPending \/ ceilingFreeSample) THEN "presentation"
    ELSE IF NormalizedCleanPass /\ handoff = "allocation" /\
            (~capacityPending \/ ceilingFreeSample) THEN "allocation"
    ELSE IF capacityEligible /\ capacityProducer THEN "capacity"
    ELSE IF retainedRefinementPending THEN "observation"
    ELSE "none"

TypeOK ==
    /\ phase \in Phases
    /\ growthPending \in BOOLEAN
    /\ coveragePending \in BOOLEAN
    /\ structuralPending \in BOOLEAN
    /\ retainedRefinementPending \in BOOLEAN
    /\ selectiveDelta \in BOOLEAN
    /\ changedCut \in BOOLEAN
    /\ rescanAfterFrame \in BOOLEAN
    /\ handoff \in Handoffs
    /\ allocationCertified \in BOOLEAN
    /\ ceilingFreeSample \in BOOLEAN
    /\ capacityPending \in BOOLEAN
    /\ capacityEligible \in BOOLEAN
    /\ capacityProducer \in BOOLEAN
    /\ capacityPopulation \in Populations
    /\ capacityCandidate \in Populations
    /\ initialCapacityPopulation \in Populations
    /\ capacityRevision \in 0..1
    /\ owner \in Owners
    /\ effect \in Effects
    /\ completedPassCount \in 0..MaximumCompletedPassCount

Init ==
    /\ phase = "completed"
    /\ growthPending \in BOOLEAN
    /\ coveragePending \in BOOLEAN
    /\ structuralPending \in BOOLEAN
    /\ retainedRefinementPending \in BOOLEAN
    /\ selectiveDelta \in BOOLEAN
    /\ changedCut \in BOOLEAN
    /\ rescanAfterFrame \in BOOLEAN
    /\ capacityEligible \in BOOLEAN
    /\ (changedCut \/ rescanAfterFrame) => capacityEligible
    /\ handoff \in Handoffs
    /\ allocationCertified \in BOOLEAN
    /\ ceilingFreeSample \in BOOLEAN
    /\ capacityPending \in BOOLEAN
    /\ ceilingFreeSample => capacityPending
    /\ capacityPending => capacityEligible
    /\ capacityProducer = capacityEligible
    /\ ceilingFreeSample =>
           handoff = "allocation" /\
           (rescanAfterFrame \/ changedCut) /\ allocationCertified
    /\ capacityPopulation \in Populations
    /\ capacityCandidate \in Populations
    /\ initialCapacityPopulation = capacityPopulation
    /\ capacityRevision = 0
    /\ owner = "none"
    /\ effect = "none"
    /\ completedPassCount = 0

SelectOwner ==
    /\ phase = "completed"
    /\ phase' = "execute"
    /\ owner' = ExpectedOwner
    /\ effect' = "none"
    /\ selectiveDelta' =
           IF ExpectedOwner \in {"allocation", "capacity"}
           THEN FALSE
           ELSE selectiveDelta
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    retainedRefinementPending,
                    changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision,
                    completedPassCount>>

ExecuteGrowth ==
    /\ phase = "execute"
    /\ owner = "growth"
    /\ phase' = "completed"
    /\ growthPending' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "drainGrowth"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<coveragePending, structuralPending, selectiveDelta, handoff,
                    allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCoverage ==
    /\ phase = "execute"
    /\ owner = "coverage"
    /\ phase' = "completed"
    /\ coveragePending' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "retryCoverage"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, structuralPending, selectiveDelta, handoff,
                    allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteStructural ==
    /\ phase = "execute"
    /\ owner = "structural"
    /\ phase' = "completed"
    /\ structuralPending' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ allocationCertified' = FALSE
    /\ ceilingFreeSample' = FALSE
    /\ capacityPending' = FALSE
    /\ capacityEligible' = FALSE
    /\ capacityProducer' = FALSE
    /\ effect' = "repairStructural"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, selectiveDelta, handoff,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecutePresentationHandoff ==
    /\ phase = "execute"
    /\ owner = "presentation"
    /\ phase' = "completed"
    /\ handoff' = "allocation"
    /\ retainedRefinementPending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "presentHandoff"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    selectiveDelta,
                    allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteAllocationHandoff ==
    /\ phase = "execute"
    /\ owner = "allocation"
    /\ phase' = "completed"
    /\ changedCut' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ selectiveDelta' = FALSE
    /\ completedPassCount' = completedPassCount + 1
    /\ capacityProducer' = capacityEligible
    /\ IF allocationCertified
          THEN /\ handoff' = "none"
               /\ allocationCertified' = TRUE
               /\ ceilingFreeSample' = FALSE
               /\ effect' = "finishHandoff"
          ELSE /\ handoff' = "allocation"
               /\ allocationCertified' = TRUE
               /\ ceilingFreeSample' = ceilingFreeSample
               /\ effect' = "allocateHandoff"
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    capacityPending,
                    capacityEligible,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCapacityCut ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ (changedCut \/ rescanAfterFrame)
    /\ phase' = "completed"
    /\ changedCut' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "presentCut"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    selectiveDelta, handoff,
                    allocationCertified, ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCapacityPopulation ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ CleanPass
    /\ capacityPopulation # capacityCandidate
    /\ phase' = "completed"
    /\ capacityPopulation' = capacityCandidate
    /\ retainedRefinementPending' = FALSE
    /\ capacityRevision' = capacityRevision + 1
    /\ effect' = "changePopulation"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    selectiveDelta, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer, capacityCandidate,
                    initialCapacityPopulation, owner>>

ExecuteCapacitySample ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ CleanPass
    /\ capacityPopulation = capacityCandidate
    /\ phase' = "completed"
    /\ capacityEligible' = FALSE
    /\ retainedRefinementPending' = FALSE
    /\ capacityPending' = FALSE
    /\ capacityProducer' = FALSE
    /\ ceilingFreeSample' = FALSE
    /\ effect' = "recordSample"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    selectiveDelta, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteRetainedObservation ==
    /\ phase = "execute"
    /\ owner = "observation"
    /\ phase' = "completed"
    /\ retainedRefinementPending' = FALSE
    /\ effect' = "retireObservation"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    selectiveDelta, changedCut, rescanAfterFrame, handoff,
                    allocationCertified, ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteNone ==
    /\ phase = "execute"
    /\ owner = "none"
    /\ phase' = "terminal"
    /\ effect' = "none"
    /\ UNCHANGED <<growthPending, coveragePending, structuralPending,
                    retainedRefinementPending,
                    selectiveDelta, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample, capacityPending,
                    capacityEligible, capacityProducer,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner,
                    completedPassCount>>

Next ==
    \/ SelectOwner
    \/ ExecuteGrowth
    \/ ExecuteCoverage
    \/ ExecuteStructural
    \/ ExecutePresentationHandoff
    \/ ExecuteAllocationHandoff
    \/ ExecuteCapacityCut
    \/ ExecuteCapacityPopulation
    \/ ExecuteCapacitySample
    \/ ExecuteRetainedObservation
    \/ ExecuteNone

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SelectOwner)
    /\ WF_vars(ExecuteGrowth)
    /\ WF_vars(ExecuteCoverage)
    /\ WF_vars(ExecuteStructural)
    /\ WF_vars(ExecutePresentationHandoff)
    /\ WF_vars(ExecuteAllocationHandoff)
    /\ WF_vars(ExecuteCapacityCut)
    /\ WF_vars(ExecuteCapacityPopulation)
    /\ WF_vars(ExecuteCapacitySample)
    /\ WF_vars(ExecuteRetainedObservation)
    /\ WF_vars(ExecuteNone)

OwnerMatchesSnapshot ==
    phase = "execute" => owner = ExpectedOwner

EffectHasExclusiveOwner ==
    /\ effect = "drainGrowth" => owner = "growth"
    /\ effect = "retryCoverage" => owner = "coverage"
    /\ effect = "repairStructural" => owner = "structural"
    /\ effect = "presentHandoff" => owner = "presentation"
    /\ effect \in {"allocateHandoff", "finishHandoff"} =>
           owner = "allocation"
    /\ effect \in {"presentCut", "changePopulation", "recordSample"} =>
           owner = "capacity"
    /\ effect = "retireObservation" => owner = "observation"

SuccessorPassStartsFresh ==
    completedPassCount > 0 /\ phase = "completed" =>
        ~changedCut /\ ~rescanAfterFrame /\ ~retainedRefinementPending

SceneWideAllocationConsumesDelta ==
    phase = "execute" /\ owner \in {"allocation", "capacity"} =>
        ~selectiveDelta

CapacityRevisionMatchesPopulation ==
    /\ capacityRevision = 0 =>
           capacityPopulation = initialCapacityPopulation
    /\ capacityRevision = 1 =>
           capacityPopulation = capacityCandidate /\
           capacityCandidate # initialCapacityPopulation

PendingCapacityHasProducer ==
    capacityPending => capacityProducer

FinishedHandoffRestoresCapacityProducer ==
    effect = "finishHandoff" /\ capacityPending => capacityProducer

StructuralFrontierPrecedesAllocation ==
    phase = "execute" /\ structuralPending /\
    ~growthPending /\ ~coveragePending => owner = "structural"

StructuralRepairInvalidatesCapacity ==
    effect = "repairStructural" =>
        /\ ~capacityPending
        /\ ~capacityEligible
        /\ ~capacityProducer
        /\ ~allocationCertified
        /\ ~retainedRefinementPending

TerminalHasNoOwner ==
    phase = "terminal" =>
        owner = "none" /\ ~structuralPending /\ ~retainedRefinementPending

CapacityRevisionNamesSemanticChange ==
    [][capacityRevision' # capacityRevision =>
        /\ owner = "capacity"
        /\ capacityPopulation' # capacityPopulation]_vars

EventuallyResolved == <> (phase = "terminal")

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => phase = "terminal"

=============================================================================
