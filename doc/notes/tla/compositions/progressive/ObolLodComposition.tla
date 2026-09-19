------------------------- MODULE ObolLodComposition -------------------------
\* Composition contract for the drawing controller's focused proof boundaries.
\*
\* The detailed admission, retained-growth, presentation-handoff, capacity,
\* structural-frontier, and point-quality models each prove a local protocol.
\* This model checks the assumptions at their shared seams: one semantic owner
\* at a time, no fallback regression while retained data is replanned, exact
\* visibility before timing or readiness, capacity-owned growth staying inside
\* its candidate, and a finite path to one terminal outcome.
\*
\* Geometry arithmetic, visual-importance scores, frame duration, and worker
\* scheduling remain outside this finite control model.

EXTENDS Naturals, TLC

CONSTANT SampleLimit

ASSUME SampleLimit > 0

Candidates == {"hero", "ordinary", "tiny"}
Significant == Candidates \ {"tiny"}
\* Source profile certification, not a workload/cardinality profile.
Profiles == {"unknown", "safe", "large"}
Modes == {"box", "aggregate", "proxy", "pop", "direct"}
Owners == {"admission", "growth", "presentation", "reconcile",
           "capacity", "structural", "point", "static", "none"}
VisibilityStates == {"unknown", "empty", "nonempty"}
Outcomes == {"none", "ready", "constrained"}
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
          owner,
          actualVisible,
          visibilityCensus,
          presentationExact,
          presentedWork,
          ceilingEffective,
          independentGrowthPending,
          capacityGrowthPending,
          capacityPending,
          structuralPending,
          pointPending,
          staticPending,
          samplesRemaining,
          allocationRevision,
          allocationCurrent,
          constraintWitness,
          outcome

vars == <<profile, capacity, mode, owner, actualVisible, visibilityCensus,
          presentationExact, presentedWork, ceilingEffective,
          independentGrowthPending, capacityGrowthPending, capacityPending,
          structuralPending, pointPending, staticPending, samplesRemaining,
          allocationRevision, allocationCurrent, constraintWitness, outcome>>

Ready == outcome = "ready"
Constrained == outcome = "constrained"
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

AllPresented == \A candidate \in Candidates: mode[candidate] # "box"

AffordableDesiredUpgrade(candidate) ==
    IF profile = "safe"
    THEN /\ mode[candidate] \in {"box", "pop"}
         /\ TotalCost(Replace(mode, candidate, "direct")) <= capacity
    ELSE /\ mode[candidate] = "box"
         /\ TotalCost(Replace(mode, candidate, "pop")) <= capacity

AnyAffordableDesiredUpgrade ==
    \E candidate \in Significant: AffordableDesiredUpgrade(candidate)

ExactVisibilityPresentation ==
    /\ presentationExact
    /\ visibilityCensus # "unknown"
    /\ (presentedWork <=> visibilityCensus = "nonempty")

PresentationSuccessor ==
    IF independentGrowthPending THEN "growth"
    ELSE IF ceilingEffective THEN "reconcile"
    ELSE IF structuralPending THEN "structural"
    ELSE IF pointPending THEN "point"
    ELSE IF capacityPending THEN "capacity"
    ELSE IF staticPending THEN "static"
    ELSE "none"

TypeOK ==
    /\ profile \in Profiles
    /\ capacity \in InitialCapacities
    /\ mode \in [Candidates -> Modes]
    /\ owner \in Owners
    /\ actualVisible \in BOOLEAN
    /\ visibilityCensus \in VisibilityStates
    /\ presentationExact \in BOOLEAN
    /\ presentedWork \in BOOLEAN
    /\ ceilingEffective \in BOOLEAN
    /\ independentGrowthPending \in BOOLEAN
    /\ capacityGrowthPending \in BOOLEAN
    /\ capacityPending \in BOOLEAN
    /\ structuralPending \in BOOLEAN
    /\ pointPending \in BOOLEAN
    /\ staticPending \in BOOLEAN
    /\ samplesRemaining \in 0..SampleLimit
    /\ allocationRevision \in 0..1
    /\ allocationCurrent \in BOOLEAN
    /\ constraintWitness \in BOOLEAN
    /\ outcome \in Outcomes

Init ==
    /\ profile = "unknown"
    /\ capacity \in InitialCapacities
    /\ mode = [candidate \in Candidates |-> "box"]
    /\ owner = "admission"
    /\ actualVisible \in BOOLEAN
    /\ visibilityCensus = "unknown"
    /\ presentationExact = FALSE
    /\ presentedWork = FALSE
    /\ ceilingEffective \in BOOLEAN
    /\ independentGrowthPending \in BOOLEAN
    /\ capacityGrowthPending \in BOOLEAN
    /\ capacityPending \in BOOLEAN
    /\ capacityGrowthPending => capacityPending
    /\ structuralPending \in BOOLEAN
    /\ pointPending \in BOOLEAN
    /\ staticPending \in BOOLEAN
    /\ samplesRemaining = SampleLimit
    /\ allocationRevision = 0
    /\ allocationCurrent = FALSE
    /\ constraintWitness = FALSE
    /\ outcome = "none"

CertifyProfile(nextProfile) ==
    /\ owner = "admission"
    /\ profile = "unknown"
    /\ nextProfile \in {"safe", "large"}
    /\ profile' = nextProfile
    /\ visibilityCensus' = IF actualVisible THEN "nonempty" ELSE "empty"
    /\ UNCHANGED <<capacity, mode, owner, actualVisible,
                    presentationExact, presentedWork, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    staticPending, samplesRemaining, allocationRevision,
                    allocationCurrent,
                    constraintWitness, outcome>>

ClassifyTiny ==
    /\ owner = "admission"
    /\ mode["tiny"] = "box"
    /\ mode' = [mode EXCEPT !["tiny"] = "aggregate"]
    /\ UNCHANGED <<profile, capacity, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

AdmitPop(candidate) ==
    LET replacement == Replace(mode, candidate, "pop") IN
    /\ owner = "admission"
    /\ candidate \in Significant
    /\ mode[candidate] = "box"
    /\ TotalCost(replacement) <= capacity
    /\ mode' = replacement
    /\ UNCHANGED <<profile, capacity, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

AdmitDirect(candidate) ==
    LET replacement == Replace(mode, candidate, "direct") IN
    /\ owner = "admission"
    /\ profile = "safe"
    /\ candidate \in Significant
    /\ mode[candidate] \in {"box", "pop"}
    /\ TotalCost(replacement) <= capacity
    /\ mode' = replacement
    /\ UNCHANGED <<profile, capacity, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

PublishPersistentProxy(candidate) ==
    LET popReplacement == Replace(mode, candidate, "pop")
        directReplacement == Replace(mode, candidate, "direct") IN
    /\ owner = "admission"
    /\ profile # "unknown"
    /\ candidate \in Significant
    /\ mode[candidate] = "box"
    /\ TotalCost(popReplacement) > capacity
    /\ (profile # "safe" \/ TotalCost(directReplacement) > capacity)
    /\ mode' = Replace(mode, candidate, "proxy")
    /\ constraintWitness' = TRUE
    /\ UNCHANGED <<profile, capacity, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent, outcome>>

CompleteAdmission ==
    /\ owner = "admission"
    /\ profile # "unknown"
    /\ AllPresented
    /\ ~AnyAffordableDesiredUpgrade
    /\ owner' = "presentation"
    /\ allocationCurrent' = TRUE
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    staticPending, samplesRemaining, allocationRevision,
                    constraintWitness, outcome>>

CompletePresentation ==
    /\ owner = "presentation"
    /\ allocationCurrent
    /\ AllPresented
    /\ visibilityCensus # "unknown"
    /\ owner' = PresentationSuccessor
    /\ presentationExact' = TRUE
    /\ presentedWork' = actualVisible
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    staticPending, samplesRemaining, allocationRevision,
                    allocationCurrent,
                    constraintWitness, outcome>>

CompleteIndependentGrowth ==
    /\ owner = "growth"
    /\ independentGrowthPending
    /\ independentGrowthPending' = FALSE
    /\ owner' = "admission"
    /\ allocationRevision' = allocationRevision + 1
    /\ allocationCurrent' = FALSE
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    constraintWitness, outcome>>

CompleteReconciliation ==
    /\ owner = "reconcile"
    /\ ceilingEffective
    /\ ceilingEffective' = FALSE
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

AbsorbCapacityGrowth ==
    /\ owner = "capacity"
    /\ capacityGrowthPending
    /\ capacityGrowthPending' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    staticPending, samplesRemaining, allocationRevision,
                    allocationCurrent,
                    constraintWitness, outcome>>

ConsumeCapacitySample ==
    /\ owner = "capacity"
    /\ ~capacityGrowthPending
    /\ capacityPending
    /\ samplesRemaining > 0
    /\ allocationCurrent
    /\ ExactVisibilityPresentation
    /\ samplesRemaining' = samplesRemaining - 1
    /\ capacityPending' = (samplesRemaining > 1)
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    structuralPending, pointPending, staticPending,
                    allocationRevision,
                    allocationCurrent, constraintWitness, outcome>>

CompleteStructuralRepair ==
    /\ owner = "structural"
    /\ structuralPending
    /\ structuralPending' = FALSE
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

CompletePointQuality ==
    /\ owner = "point"
    /\ pointPending
    /\ pointPending' = FALSE
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

CompleteStaticQuality ==
    /\ owner = "static"
    /\ staticPending
    /\ staticPending' = FALSE
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    samplesRemaining, allocationRevision, allocationCurrent,
                    constraintWitness, outcome>>

\* A rejected static-quality successor is terminal evidence only after its
\* retained framebuffer is presented exactly.  It cannot merely clear the
\* owner and let the composition publish an unwitnessed ready result.
ConstrainStaticQuality ==
    /\ owner = "static"
    /\ staticPending
    /\ staticPending' = FALSE
    /\ constraintWitness' = TRUE
    /\ owner' = "presentation"
    /\ presentationExact' = FALSE
    /\ presentedWork' = FALSE
    /\ UNCHANGED <<profile, capacity, mode, actualVisible,
                    visibilityCensus, ceilingEffective,
                    independentGrowthPending, capacityGrowthPending,
                    capacityPending, structuralPending, pointPending,
                    samplesRemaining, allocationRevision, allocationCurrent,
                    outcome>>

PublishTerminal ==
    /\ owner = "none"
    /\ outcome = "none"
    /\ allocationCurrent
    /\ ExactVisibilityPresentation
    /\ ~independentGrowthPending
    /\ ~capacityGrowthPending
    /\ ~capacityPending
    /\ ~structuralPending
    /\ ~pointPending
    /\ ~staticPending
    /\ AllPresented
    /\ outcome' = IF constraintWitness THEN "constrained" ELSE "ready"
    /\ UNCHANGED <<profile, capacity, mode, owner, actualVisible,
                    visibilityCensus, presentationExact, presentedWork,
                    ceilingEffective, independentGrowthPending,
                    capacityGrowthPending, capacityPending,
                    structuralPending, pointPending, staticPending,
                    samplesRemaining,
                    allocationRevision, allocationCurrent,
                    constraintWitness>>

Next ==
    \/ \E nextProfile \in {"safe", "large"}: CertifyProfile(nextProfile)
    \/ ClassifyTiny
    \/ \E candidate \in Significant: AdmitPop(candidate)
    \/ \E candidate \in Significant: AdmitDirect(candidate)
    \/ \E candidate \in Significant: PublishPersistentProxy(candidate)
    \/ CompleteAdmission
    \/ CompletePresentation
    \/ CompleteIndependentGrowth
    \/ CompleteReconciliation
    \/ AbsorbCapacityGrowth
    \/ ConsumeCapacitySample
    \/ CompleteStructuralRepair
    \/ CompletePointQuality
    \/ CompleteStaticQuality
    \/ ConstrainStaticQuality
    \/ PublishTerminal

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ \A nextProfile \in {"safe", "large"}:
           WF_vars(CertifyProfile(nextProfile))
    /\ WF_vars(ClassifyTiny)
    /\ \A candidate \in Significant: WF_vars(AdmitPop(candidate))
    /\ \A candidate \in Significant: WF_vars(AdmitDirect(candidate))
    /\ \A candidate \in Significant:
           WF_vars(PublishPersistentProxy(candidate))
    /\ WF_vars(CompleteAdmission)
    /\ WF_vars(CompletePresentation)
    /\ WF_vars(CompleteIndependentGrowth)
    /\ WF_vars(CompleteReconciliation)
    /\ WF_vars(AbsorbCapacityGrowth)
    /\ WF_vars(ConsumeCapacitySample)
    /\ WF_vars(CompleteStructuralRepair)
    /\ WF_vars(CompletePointQuality)
    /\ WF_vars(CompleteStaticQuality)
    /\ WF_vars(ConstrainStaticQuality)
    /\ WF_vars(PublishTerminal)

BudgetSafety == TotalCost(mode) <= capacity

VisibilityCensusIsAuthoritative ==
    /\ (visibilityCensus = "empty" => ~actualVisible)
    /\ (visibilityCensus = "nonempty" => actualVisible)

PersistentProxyHasWitness ==
    (\E candidate \in Significant: mode[candidate] = "proxy") =>
        constraintWitness

CapacitySampleRequiresExactPopulation ==
    owner = "capacity" => ExactVisibilityPresentation

TerminalIsExactAndQuiescent ==
    outcome # "none" =>
        /\ owner = "none"
        /\ ExactVisibilityPresentation
        /\ ~independentGrowthPending
        /\ ~capacityGrowthPending
        /\ ~capacityPending
        /\ ~structuralPending
        /\ ~pointPending
        /\ ~staticPending
        /\ AllPresented

TerminalHasNoStructuralBoxes ==
    outcome # "none" => \A candidate \in Candidates: mode[candidate] # "box"

ZeroWorkTerminalRequiresExactEmpty ==
    outcome # "none" /\ ~presentedWork => visibilityCensus = "empty"

CapacityGrowthHasNoSecondOwner ==
    owner = "capacity" /\ capacityGrowthPending => allocationCurrent

DirectIsAbsorbing ==
    \A candidate \in Significant:
        [](mode[candidate] = "direct" => [](mode[candidate] = "direct"))

ZeroMarginalSafePopEventuallyDirect ==
    \A candidate \in Significant:
        [](profile = "safe" /\ owner = "admission" /\
           mode[candidate] = "pop" /\
           TotalCost(Replace(mode, candidate, "direct")) = TotalCost(mode)
           => <>(mode[candidate] = "direct"))

TerminalDoesNotReopen ==
    [](outcome # "none" => [](outcome # "none"))

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyTerminal == <>(outcome # "none")

=============================================================================
