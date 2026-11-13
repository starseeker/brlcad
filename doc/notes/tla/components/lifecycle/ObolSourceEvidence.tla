----------------------- MODULE ObolSourceEvidence -----------------------
\* Exact source observation and consumption refine the source-input edge of
\* the submission/lifecycle contracts.  Several mutations and notifications
\* may precede one submission; an observed-but-unconsumed Boolean cannot
\* authenticate the next observation.  Geometry, source cardinality, endpoint
\* retirement and allocation arithmetic remain outside this component.

EXTENDS Naturals, FiniteSets, TLC

CONSTANT MaxMutations

EmptySnapshot == [population |-> 0, inventory |-> 0, visibility |-> 0]
Snapshots == [population : 0..MaxMutations,
              inventory : 0..MaxMutations, visibility : 0..MaxMutations]
NoCandidate == {}

VARIABLES source, observed, delivered, submitted, candidate, mutations, inputClosed,
          inventoryEpoch, visibilityEpoch

vars == <<source, observed, delivered, submitted, candidate, mutations, inputClosed,
          inventoryEpoch, visibilityEpoch>>

InventoryChanged(current, previous) ==
    current.population # previous.population \/
    current.inventory # previous.inventory

Init ==
    /\ source = EmptySnapshot
    /\ observed = EmptySnapshot
    /\ delivered = EmptySnapshot
    /\ submitted = EmptySnapshot
    /\ candidate = NoCandidate
    /\ mutations = 0
    /\ inputClosed = FALSE
    /\ inventoryEpoch = 0
    /\ visibilityEpoch = 0

MutateInventory ==
    /\ ~inputClosed /\ mutations < MaxMutations
    /\ source' = [source EXCEPT !.inventory = @ + 1]
    /\ mutations' = mutations + 1
    /\ UNCHANGED <<observed, delivered, submitted, candidate, inputClosed,
                    inventoryEpoch, visibilityEpoch>>

MutateVisibility ==
    /\ ~inputClosed /\ mutations < MaxMutations
    /\ source' = [source EXCEPT !.visibility = @ + 1]
    /\ mutations' = mutations + 1
    /\ UNCHANGED <<observed, delivered, submitted, candidate, inputClosed,
                    inventoryEpoch, visibilityEpoch>>

ReplacePopulation ==
    /\ ~inputClosed /\ mutations < MaxMutations
    /\ source' = [source EXCEPT !.population = @ + 1]
    /\ mutations' = mutations + 1
    /\ UNCHANGED <<observed, delivered, submitted, candidate, inputClosed,
                    inventoryEpoch, visibilityEpoch>>

PublishSnapshot(current) ==
    /\ observed' = current
    /\ inventoryEpoch' = inventoryEpoch +
        IF InventoryChanged(current, observed) THEN 1 ELSE 0
    /\ visibilityEpoch' = visibilityEpoch +
        IF current # observed /\ ~InventoryChanged(current, observed)
        THEN 1 ELSE 0

\* delivered is a ghost record of the input crossing the observation
\* interface, not another production ledger.  It makes a wrongly ignored
\* notification distinguishable from legitimate stuttering/duplicate input.
Observe ==
    /\ source # delivered
    /\ delivered' = source
    /\ PublishSnapshot(source)
    /\ UNCHANGED <<source, submitted, candidate, mutations, inputClosed>>

\* Submission also observes: it cannot rely on delivery of a separate host
\* notification.  Publishing an already observed snapshot is an inert step.
BeginSubmission ==
    /\ candidate = NoCandidate
    /\ submitted # source
    /\ delivered' = source
    /\ PublishSnapshot(source)
    /\ candidate' = {source}
    /\ UNCHANGED <<source, submitted, mutations, inputClosed>>

FinishSubmission ==
    /\ candidate # NoCandidate
    /\ submitted' = IF observed \in candidate THEN observed ELSE submitted
    /\ candidate' = NoCandidate
    /\ UNCHANGED <<source, observed, delivered, mutations, inputClosed,
                    inventoryEpoch, visibilityEpoch>>

CloseInputs ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<source, observed, delivered, submitted, candidate, mutations,
                    inventoryEpoch, visibilityEpoch>>

Next == MutateInventory \/ MutateVisibility \/ ReplacePopulation \/
        Observe \/ BeginSubmission \/ FinishSubmission \/ CloseInputs

Spec == Init /\ [][Next]_vars /\ WF_vars(Observe) /\
        WF_vars(BeginSubmission) /\ WF_vars(FinishSubmission)

TypeOK ==
    /\ source \in Snapshots /\ observed \in Snapshots /\ submitted \in Snapshots
    /\ delivered \in Snapshots
    /\ candidate \subseteq Snapshots /\ Cardinality(candidate) <= 1
    /\ mutations \in 0..MaxMutations /\ inputClosed \in BOOLEAN
    /\ inventoryEpoch \in 0..MaxMutations
    /\ visibilityEpoch \in 0..MaxMutations

\* These action properties reject an unkeyed notification latch, publication
\* from the old submitted baseline, and a stale submission consuming the
\* newest observed snapshot.  Unchanged notifications are stuttering steps.
ObservedInputIsDelivered == observed = delivered

ObservationPublishesCurrentInput ==
    [][Observe => observed' = source]_vars

ExactDomainPublication ==
    [][(observed' # observed) =>
        IF InventoryChanged(observed', observed)
        THEN inventoryEpoch' = inventoryEpoch + 1 /\
             visibilityEpoch' = visibilityEpoch
        ELSE visibilityEpoch' = visibilityEpoch + 1 /\
             inventoryEpoch' = inventoryEpoch]_vars

UnchangedEvidenceIsInert ==
    [][(observed' = observed) => inventoryEpoch' = inventoryEpoch /\
                                  visibilityEpoch' = visibilityEpoch]_vars

ConsumptionMatchesObservation ==
    [][(submitted' # submitted) => observed \in candidate /\
                                   submitted' = observed]_vars

Terminal == inputClosed /\ candidate = NoCandidate /\
            source = observed /\ source = submitted
DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal
EventuallyConsumesClosedInput == inputClosed ~> Terminal

=======================================================================
