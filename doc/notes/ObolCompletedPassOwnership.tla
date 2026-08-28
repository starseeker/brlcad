----------------------- MODULE ObolCompletedPassOwnership ----------------------
\* Multi-pass successor transaction for a completed progressive-planning pass.
\*
\* Selection is made from one immutable completion snapshot before any owner
\* mutates evidence.  Resident growth and coverage have precedence over
\* presentation handoff and generic capacity calibration.  Every effect which
\* starts a successor pass clears the complete pass-annotation transaction.
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
Owners == {"none", "growth", "coverage", "presentation", "allocation",
           "capacity"}
Effects == {"none", "drainGrowth", "retryCoverage", "presentHandoff",
            "allocateHandoff", "finishHandoff", "presentCut",
            "changePopulation", "recordSample"}
Populations == 0..1
\* At most growth, coverage, three handoff effects (presentation, allocation,
\* finish), and three capacity effects can precede terminal selection in this
\* symmetry-reduced transaction.
MaximumCompletedPassCount == 8

VARIABLES phase,
          growthPending,
          coveragePending,
          changedCut,
          rescanAfterFrame,
          handoff,
          allocationCertified,
          ceilingFreeSample,
          capacityEligible,
          capacityPopulation,
          capacityCandidate,
          initialCapacityPopulation,
          capacityRevision,
          owner,
          effect,
          completedPassCount

vars == <<phase, growthPending, coveragePending, changedCut,
          rescanAfterFrame, handoff, allocationCertified, ceilingFreeSample,
          capacityEligible,
          capacityPopulation, capacityCandidate, initialCapacityPopulation,
          capacityRevision, owner, effect, completedPassCount>>

CleanPass == ~changedCut /\ ~rescanAfterFrame
NormalizedCleanPass == CleanPass \/ ceilingFreeSample

ExpectedOwner ==
    IF growthPending THEN "growth"
    ELSE IF coveragePending THEN "coverage"
    ELSE IF NormalizedCleanPass /\ handoff = "presentation" THEN "presentation"
    ELSE IF NormalizedCleanPass /\ handoff = "allocation" THEN "allocation"
    ELSE IF capacityEligible THEN "capacity"
    ELSE "none"

TypeOK ==
    /\ phase \in Phases
    /\ growthPending \in BOOLEAN
    /\ coveragePending \in BOOLEAN
    /\ changedCut \in BOOLEAN
    /\ rescanAfterFrame \in BOOLEAN
    /\ handoff \in Handoffs
    /\ allocationCertified \in BOOLEAN
    /\ ceilingFreeSample \in BOOLEAN
    /\ capacityEligible \in BOOLEAN
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
    /\ changedCut \in BOOLEAN
    /\ rescanAfterFrame \in BOOLEAN
    /\ capacityEligible \in BOOLEAN
    /\ (changedCut \/ rescanAfterFrame) => capacityEligible
    /\ handoff \in Handoffs
    /\ allocationCertified \in BOOLEAN
    /\ ceilingFreeSample \in BOOLEAN
    /\ ceilingFreeSample =>
           handoff = "allocation" /\ rescanAfterFrame /\ allocationCertified
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
    /\ UNCHANGED <<growthPending, coveragePending, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision,
                    completedPassCount>>

ExecuteGrowth ==
    /\ phase = "execute"
    /\ owner = "growth"
    /\ phase' = "completed"
    /\ growthPending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "drainGrowth"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<coveragePending, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCoverage ==
    /\ phase = "execute"
    /\ owner = "coverage"
    /\ phase' = "completed"
    /\ coveragePending' = FALSE
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "retryCoverage"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecutePresentationHandoff ==
    /\ phase = "execute"
    /\ owner = "presentation"
    /\ phase' = "completed"
    /\ handoff' = "allocation"
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "presentHandoff"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteAllocationHandoff ==
    /\ phase = "execute"
    /\ owner = "allocation"
    /\ phase' = "completed"
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ completedPassCount' = completedPassCount + 1
    /\ IF allocationCertified
          THEN /\ handoff' = "none"
               /\ allocationCertified' = TRUE
               /\ effect' = "finishHandoff"
          ELSE /\ handoff' = "allocation"
               /\ allocationCertified' = TRUE
               /\ effect' = "allocateHandoff"
    /\ UNCHANGED <<growthPending, coveragePending, ceilingFreeSample,
                    capacityEligible,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCapacityCut ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ (changedCut \/ rescanAfterFrame)
    /\ phase' = "completed"
    /\ changedCut' = FALSE
    /\ rescanAfterFrame' = FALSE
    /\ effect' = "presentCut"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, handoff,
                    allocationCertified, ceilingFreeSample, capacityEligible,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteCapacityPopulation ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ CleanPass
    /\ capacityPopulation # capacityCandidate
    /\ phase' = "completed"
    /\ capacityPopulation' = capacityCandidate
    /\ capacityRevision' = capacityRevision + 1
    /\ effect' = "changePopulation"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityCandidate,
                    initialCapacityPopulation, owner>>

ExecuteCapacitySample ==
    /\ phase = "execute"
    /\ owner = "capacity"
    /\ CleanPass
    /\ capacityPopulation = capacityCandidate
    /\ phase' = "completed"
    /\ capacityEligible' = FALSE
    /\ effect' = "recordSample"
    /\ completedPassCount' = completedPassCount + 1
    /\ UNCHANGED <<growthPending, coveragePending, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner>>

ExecuteNone ==
    /\ phase = "execute"
    /\ owner = "none"
    /\ phase' = "terminal"
    /\ effect' = "none"
    /\ UNCHANGED <<growthPending, coveragePending, changedCut,
                    rescanAfterFrame, handoff, allocationCertified,
                    ceilingFreeSample,
                    capacityEligible, capacityPopulation, capacityCandidate,
                    initialCapacityPopulation, capacityRevision, owner,
                    completedPassCount>>

Next ==
    \/ SelectOwner
    \/ ExecuteGrowth
    \/ ExecuteCoverage
    \/ ExecutePresentationHandoff
    \/ ExecuteAllocationHandoff
    \/ ExecuteCapacityCut
    \/ ExecuteCapacityPopulation
    \/ ExecuteCapacitySample
    \/ ExecuteNone

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SelectOwner)
    /\ WF_vars(ExecuteGrowth)
    /\ WF_vars(ExecuteCoverage)
    /\ WF_vars(ExecutePresentationHandoff)
    /\ WF_vars(ExecuteAllocationHandoff)
    /\ WF_vars(ExecuteCapacityCut)
    /\ WF_vars(ExecuteCapacityPopulation)
    /\ WF_vars(ExecuteCapacitySample)
    /\ WF_vars(ExecuteNone)

OwnerMatchesSnapshot ==
    phase = "execute" => owner = ExpectedOwner

EffectHasExclusiveOwner ==
    /\ effect = "drainGrowth" => owner = "growth"
    /\ effect = "retryCoverage" => owner = "coverage"
    /\ effect = "presentHandoff" => owner = "presentation"
    /\ effect \in {"allocateHandoff", "finishHandoff"} =>
           owner = "allocation"
    /\ effect \in {"presentCut", "changePopulation", "recordSample"} =>
           owner = "capacity"

SuccessorPassStartsFresh ==
    completedPassCount > 0 /\ phase = "completed" =>
        ~changedCut /\ ~rescanAfterFrame

CapacityRevisionMatchesPopulation ==
    /\ capacityRevision = 0 =>
           capacityPopulation = initialCapacityPopulation
    /\ capacityRevision = 1 =>
           capacityPopulation = capacityCandidate /\
           capacityCandidate # initialCapacityPopulation

TerminalHasNoOwner ==
    phase = "terminal" => owner = "none"

CapacityRevisionNamesSemanticChange ==
    [][capacityRevision' # capacityRevision =>
        /\ owner = "capacity"
        /\ capacityPopulation' # capacityPopulation]_vars

EventuallyResolved == <> (phase = "terminal")

=============================================================================
