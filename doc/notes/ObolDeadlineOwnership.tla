------------------------- MODULE ObolDeadlineOwnership -------------------------
\* Deadline-successor ownership below the canonical progressive pipeline.
\*
\* A failed presentation has exactly one successor.  Strict progress on the
\* same finite presentation transaction permits a retry.  Otherwise an active
\* population producer retains its cursor and eventually supplies the next
\* frame.  Capacity recovery is legal only after both owners are quiescent.
\* The model deliberately omits render cost and geometry quality; those are
\* numeric and graphical proof obligations.

EXTENDS Naturals, TLC

CONSTANT MaxPopulationWork, MaxTransactionWork

Successors == {"none", "transaction", "population", "capacity"}

VARIABLES populationRemaining,
          transactionRemaining,
          renderPending,
          capacityRecoveryPending,
          capacityRevision,
          lastSuccessor,
          priorPopulationRemaining,
          priorTransactionRemaining,
          priorCapacityRevision

vars == <<populationRemaining, transactionRemaining, renderPending,
          capacityRecoveryPending, capacityRevision, lastSuccessor,
          priorPopulationRemaining, priorTransactionRemaining,
          priorCapacityRevision>>

TypeOK ==
    /\ populationRemaining \in 0..MaxPopulationWork
    /\ transactionRemaining \in 0..MaxTransactionWork
    /\ renderPending \in BOOLEAN
    /\ capacityRecoveryPending \in BOOLEAN
    /\ capacityRevision \in 0..1
    /\ lastSuccessor \in Successors
    /\ priorPopulationRemaining \in 0..MaxPopulationWork
    /\ priorTransactionRemaining \in 0..MaxTransactionWork
    /\ priorCapacityRevision \in 0..1

Init ==
    /\ populationRemaining \in 0..MaxPopulationWork
    /\ transactionRemaining \in 0..MaxTransactionWork
    /\ renderPending = (transactionRemaining > 0)
    /\ capacityRecoveryPending = FALSE
    /\ capacityRevision = 0
    /\ lastSuccessor = "none"
    /\ priorPopulationRemaining = populationRemaining
    /\ priorTransactionRemaining = transactionRemaining
    /\ priorCapacityRevision = capacityRevision

PresentationDeadline ==
    /\ renderPending
    /\ priorPopulationRemaining' = populationRemaining
    /\ priorTransactionRemaining' = transactionRemaining
    /\ priorCapacityRevision' = capacityRevision
    /\ IF transactionRemaining > 0
          THEN /\ transactionRemaining' = transactionRemaining - 1
               /\ renderPending' = (transactionRemaining > 1)
               /\ lastSuccessor' = "transaction"
               /\ UNCHANGED <<populationRemaining,
                               capacityRecoveryPending, capacityRevision>>
          ELSE IF populationRemaining > 0
          THEN /\ renderPending' = FALSE
               /\ lastSuccessor' = "population"
               /\ UNCHANGED <<populationRemaining, transactionRemaining,
                               capacityRecoveryPending, capacityRevision>>
          ELSE /\ ~capacityRecoveryPending
               /\ capacityRevision = 0
               /\ capacityRecoveryPending' = TRUE
               /\ capacityRevision' = 1
               /\ renderPending' = FALSE
               /\ lastSuccessor' = "capacity"
               /\ UNCHANGED <<populationRemaining, transactionRemaining>>

AdvancePopulation ==
    /\ populationRemaining > 0
    /\ ~renderPending
    /\ ~capacityRecoveryPending
    /\ populationRemaining' = populationRemaining - 1
    /\ renderPending' = (populationRemaining = 1)
    /\ lastSuccessor' = "none"
    /\ UNCHANGED <<transactionRemaining, capacityRecoveryPending,
                    capacityRevision, priorPopulationRemaining,
                    priorTransactionRemaining, priorCapacityRevision>>

CompleteCapacityRecovery ==
    /\ capacityRecoveryPending
    /\ capacityRecoveryPending' = FALSE
    /\ lastSuccessor' = "none"
    /\ UNCHANGED <<populationRemaining, transactionRemaining, renderPending,
                    capacityRevision, priorPopulationRemaining,
                    priorTransactionRemaining, priorCapacityRevision>>

Next ==
    \/ PresentationDeadline
    \/ AdvancePopulation
    \/ CompleteCapacityRecovery

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PresentationDeadline)
    /\ WF_vars(AdvancePopulation)
    /\ WF_vars(CompleteCapacityRecovery)

PopulationDeadlinePreservesCursor ==
    lastSuccessor = "population" =>
        /\ populationRemaining = priorPopulationRemaining
        /\ transactionRemaining = priorTransactionRemaining
        /\ capacityRevision = priorCapacityRevision
        /\ ~renderPending

TransactionRetryMakesStrictProgress ==
    lastSuccessor = "transaction" =>
        /\ transactionRemaining < priorTransactionRemaining
        /\ populationRemaining = priorPopulationRemaining
        /\ capacityRevision = priorCapacityRevision

CapacityRecoveryRequiresQuiescence ==
    lastSuccessor = "capacity" =>
        populationRemaining = 0 /\ transactionRemaining = 0

ActiveStateHasOwner ==
    renderPending \/ populationRemaining > 0 \/ capacityRecoveryPending =>
        \/ ENABLED PresentationDeadline
        \/ ENABLED AdvancePopulation
        \/ ENABLED CompleteCapacityRecovery

EventuallyTerminal ==
    <> (~renderPending /\ populationRemaining = 0 /\
        transactionRemaining = 0 /\ ~capacityRecoveryPending)

=============================================================================
