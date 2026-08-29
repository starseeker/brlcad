------------------------- MODULE ObolDeadlineOwnership -------------------------
\* Deadline-successor ownership below the canonical progressive pipeline.
\*
\* A failed presentation has exactly one successor.  Strict progress on the
\* same finite presentation transaction permits a retry.  Otherwise an active
\* population producer retains its cursor and eventually supplies the next
\* frame.  Once both are quiescent, an already active bounded capacity search
\* owns the miss; generic renderer-ceiling recovery is legal only when no such
\* search exists.  Selecting both successors was the 150k OSMesa livelock: the
\* generic handoff changed the population which the bounded search was trying
\* to classify, so each controller restarted the other.
\*
\* Search work is deliberately abstract.  Each unit stands for one strictly
\* narrower candidate.  ObolCapacitySearch proves the numeric bracket and
\* sample bound; this model proves cross-controller ownership.

EXTENDS Naturals, TLC

CONSTANT MaxPopulationWork, MaxTransactionWork, MaxSearchWork

Successors == {"none", "transaction", "population", "search", "generic"}

VARIABLES populationRemaining,
          transactionRemaining,
          searchRemaining,
          renderPending,
          searchRecoveryPending,
          genericRecoveryPending,
          capacityRevision,
          lastSuccessor,
          priorPopulationRemaining,
          priorTransactionRemaining,
          priorSearchRemaining,
          priorCapacityRevision

vars == <<populationRemaining, transactionRemaining, searchRemaining,
          renderPending, searchRecoveryPending, genericRecoveryPending,
          capacityRevision, lastSuccessor, priorPopulationRemaining,
          priorTransactionRemaining, priorSearchRemaining,
          priorCapacityRevision>>

TypeOK ==
    /\ populationRemaining \in 0..MaxPopulationWork
    /\ transactionRemaining \in 0..MaxTransactionWork
    /\ searchRemaining \in 0..MaxSearchWork
    /\ renderPending \in BOOLEAN
    /\ searchRecoveryPending \in BOOLEAN
    /\ genericRecoveryPending \in BOOLEAN
    /\ capacityRevision \in 0..(MaxSearchWork + 1)
    /\ lastSuccessor \in Successors
    /\ priorPopulationRemaining \in 0..MaxPopulationWork
    /\ priorTransactionRemaining \in 0..MaxTransactionWork
    /\ priorSearchRemaining \in 0..MaxSearchWork
    /\ priorCapacityRevision \in 0..(MaxSearchWork + 1)

Init ==
    /\ populationRemaining \in 0..MaxPopulationWork
    /\ transactionRemaining \in 0..MaxTransactionWork
    /\ searchRemaining \in 0..MaxSearchWork
    /\ renderPending = (transactionRemaining > 0 \/
          (transactionRemaining = 0 /\ populationRemaining = 0 /\
           searchRemaining > 0))
    /\ searchRecoveryPending = FALSE
    /\ genericRecoveryPending = FALSE
    /\ capacityRevision = 0
    /\ lastSuccessor = "none"
    /\ priorPopulationRemaining = populationRemaining
    /\ priorTransactionRemaining = transactionRemaining
    /\ priorSearchRemaining = searchRemaining
    /\ priorCapacityRevision = capacityRevision

PresentationDeadline ==
    /\ renderPending
    /\ ~searchRecoveryPending
    /\ ~genericRecoveryPending
    /\ priorPopulationRemaining' = populationRemaining
    /\ priorTransactionRemaining' = transactionRemaining
    /\ priorSearchRemaining' = searchRemaining
    /\ priorCapacityRevision' = capacityRevision
    /\ IF transactionRemaining > 0
          THEN /\ transactionRemaining' = transactionRemaining - 1
               /\ renderPending' = (transactionRemaining > 1 \/
                    (transactionRemaining = 1 /\
                     populationRemaining = 0 /\ searchRemaining > 0))
               /\ lastSuccessor' = "transaction"
               /\ UNCHANGED <<populationRemaining, searchRemaining,
                               searchRecoveryPending,
                               genericRecoveryPending, capacityRevision>>
          ELSE IF populationRemaining > 0
          THEN /\ renderPending' = FALSE
               /\ lastSuccessor' = "population"
               /\ UNCHANGED <<populationRemaining, transactionRemaining,
                               searchRemaining, searchRecoveryPending,
                               genericRecoveryPending, capacityRevision>>
          ELSE IF searchRemaining > 0
          THEN /\ searchRecoveryPending' = TRUE
               /\ genericRecoveryPending' = FALSE
               /\ renderPending' = FALSE
               /\ lastSuccessor' = "search"
               /\ UNCHANGED <<populationRemaining, transactionRemaining,
                               searchRemaining, capacityRevision>>
          ELSE /\ searchRecoveryPending' = FALSE
               /\ genericRecoveryPending' = TRUE
               /\ renderPending' = FALSE
               /\ lastSuccessor' = "generic"
               /\ UNCHANGED <<populationRemaining, transactionRemaining,
                               searchRemaining, capacityRevision>>

AdvancePopulation ==
    /\ populationRemaining > 0
    /\ ~renderPending
    /\ ~searchRecoveryPending
    /\ ~genericRecoveryPending
    /\ populationRemaining' = populationRemaining - 1
    /\ renderPending' = (populationRemaining = 1)
    /\ lastSuccessor' = "none"
    /\ UNCHANGED <<transactionRemaining, searchRemaining,
                    searchRecoveryPending, genericRecoveryPending,
                    capacityRevision, priorPopulationRemaining,
                    priorTransactionRemaining, priorSearchRemaining,
                    priorCapacityRevision>>

AdvanceCapacitySearch ==
    /\ searchRecoveryPending
    /\ searchRemaining > 0
    /\ searchRecoveryPending' = FALSE
    /\ genericRecoveryPending' = FALSE
    /\ searchRemaining' = searchRemaining - 1
    /\ capacityRevision' = capacityRevision + 1
    /\ renderPending' = (searchRemaining > 1)
    /\ lastSuccessor' = "none"
    /\ UNCHANGED <<populationRemaining, transactionRemaining,
                    priorPopulationRemaining, priorTransactionRemaining,
                    priorSearchRemaining, priorCapacityRevision>>

CompleteGenericRecovery ==
    /\ genericRecoveryPending
    /\ genericRecoveryPending' = FALSE
    /\ capacityRevision' = capacityRevision + 1
    /\ lastSuccessor' = "none"
    /\ UNCHANGED <<populationRemaining, transactionRemaining,
                    searchRemaining, renderPending, searchRecoveryPending,
                    priorPopulationRemaining, priorTransactionRemaining,
                    priorSearchRemaining, priorCapacityRevision>>

Next ==
    \/ PresentationDeadline
    \/ AdvancePopulation
    \/ AdvanceCapacitySearch
    \/ CompleteGenericRecovery

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PresentationDeadline)
    /\ WF_vars(AdvancePopulation)
    /\ WF_vars(AdvanceCapacitySearch)
    /\ WF_vars(CompleteGenericRecovery)

RecoveryOwnersAreExclusive ==
    ~(searchRecoveryPending /\ genericRecoveryPending)

PopulationDeadlinePreservesCursor ==
    lastSuccessor = "population" =>
        /\ populationRemaining = priorPopulationRemaining
        /\ transactionRemaining = priorTransactionRemaining
        /\ searchRemaining = priorSearchRemaining
        /\ capacityRevision = priorCapacityRevision
        /\ ~renderPending

TransactionRetryMakesStrictProgress ==
    lastSuccessor = "transaction" =>
        /\ transactionRemaining < priorTransactionRemaining
        /\ populationRemaining = priorPopulationRemaining
        /\ searchRemaining = priorSearchRemaining
        /\ capacityRevision = priorCapacityRevision

SearchOwnsDeadline ==
    lastSuccessor = "search" =>
        /\ searchRecoveryPending
        /\ ~genericRecoveryPending
        /\ populationRemaining = 0
        /\ transactionRemaining = 0
        /\ searchRemaining > 0

GenericRequiresNoSearch ==
    lastSuccessor = "generic" =>
        /\ genericRecoveryPending
        /\ ~searchRecoveryPending
        /\ populationRemaining = 0
        /\ transactionRemaining = 0
        /\ searchRemaining = 0

ActiveStateHasOwner ==
    renderPending \/ populationRemaining > 0 \/ searchRecoveryPending \/
        genericRecoveryPending =>
        \/ ENABLED PresentationDeadline
        \/ ENABLED AdvancePopulation
        \/ ENABLED AdvanceCapacitySearch
        \/ ENABLED CompleteGenericRecovery

EventuallyTerminal ==
    <> (~renderPending /\ populationRemaining = 0 /\
        transactionRemaining = 0 /\ searchRemaining = 0 /\
        ~searchRecoveryPending /\ ~genericRecoveryPending)

=============================================================================
