------------------------ MODULE ObolLodColdPreview -------------------------
\* Bounded control-plane model for a cold large-mesh presentation.  It models
\* representation ownership, not triangle topology or wall-clock time.
\*
\* A source-order page may arrive quickly but is not globally representative.
\* Such a sample is diagnostic only: it may never replace the standing
\* overview.  Only a coverage-certified preview may do that.  A producer is
\* admitted only when its bounded working set fits remaining headroom; a
\* refusal leaves the overview visible and reports an explicit constraint.
\*
\* Binding is occurrence-local.  A coverage preview may use the adaptive mesh
\* result channel, but it is not a resident progressive binding.  Readiness
\* requires every occurrence to bind the completed hierarchy, including a
\* provisional occurrence which happened to publish the shared asset first.

EXTENDS Naturals, TLC

CONSTANT MaxEpoch

Candidates == {"prominent", "ordinary"}
Overview == "overview"
Coverage == "coverage"
Mesh == "mesh"
Aggregate == "aggregate"
Presentations == {Overview, Coverage, Mesh, Aggregate}
BindingStates == {"coverage", "resident"}
GenerationKinds == {"whole", "spatial"}

VARIABLES epoch, inputOpen, headroom, producerWorkingSet, producerPending,
          sampleReady, coverageCertified, generationKind, preview, binding,
          framePending, constrained

vars == <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
          sampleReady, coverageCertified, generationKind, preview, binding,
          framePending, constrained>>

TypeOK ==
    /\ epoch \in 0..MaxEpoch
    /\ inputOpen \in BOOLEAN
    /\ headroom \in 0..3
    /\ producerWorkingSet \in 1..3
    /\ producerPending \in BOOLEAN
    /\ sampleReady \in BOOLEAN
    /\ coverageCertified \in BOOLEAN
    /\ generationKind \in GenerationKinds
    /\ preview \in Presentations
    /\ binding \in [Candidates -> BindingStates]
    /\ framePending \in BOOLEAN
    /\ constrained \in BOOLEAN

Init ==
    /\ epoch = 0
    /\ inputOpen = FALSE
    /\ headroom \in 0..3
    /\ producerWorkingSet \in 1..3
    /\ producerPending = FALSE
    /\ sampleReady = FALSE
    /\ coverageCertified = FALSE
    /\ generationKind = "whole"
    /\ preview = Overview
    /\ binding = [candidate \in Candidates |-> "coverage"]
    /\ framePending = FALSE
    /\ constrained = FALSE

StartProducer ==
    /\ ~inputOpen /\ ~producerPending /\ ~constrained
    /\ ~coverageCertified /\ ~framePending
    /\ producerWorkingSet <= headroom
    /\ producerPending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, sampleReady,
                    coverageCertified, generationKind, preview, binding,
                    framePending, constrained>>

RefuseProducer ==
    /\ ~inputOpen /\ ~producerPending /\ ~coverageCertified
    /\ producerWorkingSet > headroom
    /\ constrained' = TRUE
    /\ preview' = Aggregate
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
                    sampleReady, coverageCertified, generationKind, binding>>

PublishSample ==
    /\ producerPending /\ ~sampleReady
    /\ sampleReady' = TRUE
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
                    coverageCertified, generationKind, preview, binding,
                    constrained>>

CertifyCoverage ==
    /\ producerPending /\ sampleReady /\ ~coverageCertified
    /\ coverageCertified' = TRUE
    /\ preview' = Coverage
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
                    sampleReady, generationKind, binding, constrained>>

CompleteSpatialHierarchy ==
    /\ producerPending /\ coverageCertified
    /\ generationKind = "whole"
    /\ generationKind' = "spatial"
    /\ producerPending' = FALSE
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, sampleReady,
                    coverageCertified, preview, binding, constrained>>

PublishMesh(candidate) ==
    /\ candidate \in Candidates
    /\ coverageCertified /\ generationKind = "spatial"
    /\ binding[candidate] = "coverage"
    /\ binding' = [binding EXCEPT ![candidate] = "resident"]
    /\ preview' = Mesh
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
                    sampleReady, coverageCertified, generationKind, constrained>>

Constrain ==
    /\ ~inputOpen /\ ~constrained
    /\ (coverageCertified \/ producerWorkingSet > headroom)
    /\ constrained' = TRUE
    /\ producerPending' = FALSE
    /\ preview' = Aggregate
    /\ framePending' = TRUE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, sampleReady,
                    coverageCertified, generationKind, binding>>

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ UNCHANGED <<epoch, inputOpen, headroom, producerWorkingSet, producerPending,
                    sampleReady, coverageCertified, generationKind, preview,
                    binding, constrained>>

BeginInput ==
    /\ ~inputOpen /\ ~framePending /\ epoch < MaxEpoch
    /\ epoch' = epoch + 1
    /\ inputOpen' = TRUE
    /\ UNCHANGED <<headroom, producerWorkingSet, producerPending, sampleReady,
                    coverageCertified, generationKind, preview, binding,
                    framePending, constrained>>

EndInput ==
    /\ inputOpen
    /\ inputOpen' = FALSE
    /\ UNCHANGED <<epoch, headroom, producerWorkingSet, producerPending, sampleReady,
                    coverageCertified, generationKind, preview, binding,
                    framePending, constrained>>

Next ==
    \/ StartProducer
    \/ RefuseProducer
    \/ PublishSample
    \/ CertifyCoverage
    \/ CompleteSpatialHierarchy
    \/ \E candidate \in Candidates: PublishMesh(candidate)
    \/ Constrain
    \/ CompleteFrame
    \/ BeginInput
    \/ EndInput

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(StartProducer)
    /\ WF_vars(RefuseProducer)
    /\ WF_vars(PublishSample)
    /\ WF_vars(CertifyCoverage)
    /\ WF_vars(CompleteSpatialHierarchy)
    /\ \A candidate \in Candidates: WF_vars(PublishMesh(candidate))
    /\ WF_vars(Constrain)
    /\ WF_vars(CompleteFrame)

\* An unrepresentative prefix never becomes the complete object's only visual.
\* The only pre-certificate replacement is the explicit resource-constrained
\* aggregate, which is complete coverage by construction rather than a sample.
UncertifiedSampleRetainsOverview ==
    ~coverageCertified /\ ~constrained => preview = Overview

WorkingSetNeverOvercommitted == producerPending => producerWorkingSet <= headroom

\* Every state retains an immediately displayable representation.
AlwaysDisplayable == preview \in Presentations

ResidentBindingRequiresSpatial ==
    (\E c \in Candidates: binding[c] = "resident")
    => generationKind = "spatial"

\* A finite input sequence has no quiet ownerless middle state.
QuietIncompleteHasOwner ==
    ~inputOpen /\ ~framePending /\ ~constrained /\ ~coverageCertified
    => producerPending \/ ENABLED StartProducer \/ ENABLED RefuseProducer

Stable ==
    /\ epoch = MaxEpoch
    /\ ~inputOpen /\ ~producerPending /\ ~framePending
    /\ (constrained \/
        (coverageCertified /\ generationKind = "spatial" /\
         \A c \in Candidates: binding[c] = "resident"))

StableNeverUsesOverview == Stable => preview \in {Coverage, Mesh, Aggregate}

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Stable

EventuallyStableAfterFinalInput ==
    []((epoch = MaxEpoch /\ ~inputOpen) => <>Stable)

=============================================================================
