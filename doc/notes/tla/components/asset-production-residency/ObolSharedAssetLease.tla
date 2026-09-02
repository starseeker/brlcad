------------------------- MODULE ObolSharedAssetLease ------------------------
\* Multi-consumer lifetime contract for one immutable shared asset build.
\*
\* Per-view demand is not asset-production lifetime.  Either view may close
\* while a shared build or result exists; cancellation is legal only after no
\* live consumer remains.  A surviving consumer retains the producer/result
\* witness and can still receive the completed asset.  The two-consumer bound
\* covers the ownership race without modeling geometry or reference counts.

EXTENDS FiniteSets, Naturals, TLC

Consumers == {"left", "right"}
ProducerStates == {"idle", "building", "result", "cancelled"}

CONSTANT MaxConsumerRejoins

VARIABLES consumerOpen,
          published,
          rejoinCount,
          producerState

vars == <<consumerOpen, published, rejoinCount, producerState>>

LiveConsumers == {consumer \in Consumers: consumerOpen[consumer]}
UnsatisfiedConsumers ==
    {consumer \in Consumers:
        consumerOpen[consumer] /\ ~published[consumer]}

TypeOK ==
    /\ consumerOpen \in [Consumers -> BOOLEAN]
    /\ published \in [Consumers -> BOOLEAN]
    /\ rejoinCount \in [Consumers -> 0..MaxConsumerRejoins]
    /\ producerState \in ProducerStates

Init ==
    /\ consumerOpen = [consumer \in Consumers |-> TRUE]
    /\ published = [consumer \in Consumers |-> FALSE]
    /\ rejoinCount = [consumer \in Consumers |-> 0]
    /\ producerState = "idle"

SubmitSharedBuild ==
    /\ producerState = "idle"
    /\ UnsatisfiedConsumers # {}
    /\ producerState' = "building"
    /\ UNCHANGED <<consumerOpen, published, rejoinCount>>

CompleteSharedBuild ==
    /\ producerState = "building"
    /\ producerState' = "result"
    /\ UNCHANGED <<consumerOpen, published, rejoinCount>>

PublishTo(consumer) ==
    /\ producerState = "result"
    /\ consumer \in UnsatisfiedConsumers
    /\ published' = [published EXCEPT ![consumer] = TRUE]
    /\ UNCHANGED <<consumerOpen, rejoinCount, producerState>>

RetireSharedResult ==
    /\ producerState = "result"
    /\ UnsatisfiedConsumers = {}
    /\ producerState' = "idle"
    /\ UNCHANGED <<consumerOpen, published, rejoinCount>>

\* Session closure is an environment action.  It retires only that view's
\* lease; it cannot mutate another consumer's evidence or the shared build.
CloseConsumer(consumer) ==
    /\ consumer \in Consumers
    /\ consumerOpen[consumer]
    /\ consumerOpen' = [consumerOpen EXCEPT ![consumer] = FALSE]
    /\ published' = [published EXCEPT ![consumer] = FALSE]
    /\ UNCHANGED <<rejoinCount, producerState>>

\* A consumer slot may represent a later view lifetime.  Rejoining is bounded
\* only to keep TLC's environment finite; it can happen during a build, while
\* a result waits for consumers, or after the preceding unowned build was
\* cancelled.  A new lifetime never inherits the old lifetime's publication.
JoinConsumer(consumer) ==
    /\ consumer \in Consumers
    /\ ~consumerOpen[consumer]
    /\ rejoinCount[consumer] < MaxConsumerRejoins
    /\ consumerOpen' = [consumerOpen EXCEPT ![consumer] = TRUE]
    /\ published' = [published EXCEPT ![consumer] = FALSE]
    /\ rejoinCount' =
          [rejoinCount EXCEPT ![consumer] = @ + 1]
    /\ producerState' =
          IF producerState = "cancelled" THEN "idle" ELSE producerState

CancelUnownedBuild ==
    /\ producerState = "building"
    /\ LiveConsumers = {}
    /\ producerState' = "cancelled"
    /\ UNCHANGED <<consumerOpen, published, rejoinCount>>

Next ==
    \/ SubmitSharedBuild
    \/ CompleteSharedBuild
    \/ \E consumer \in Consumers: PublishTo(consumer)
    \/ RetireSharedResult
    \/ \E consumer \in Consumers: CloseConsumer(consumer)
    \/ \E consumer \in Consumers: JoinConsumer(consumer)
    \/ CancelUnownedBuild

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SubmitSharedBuild)
    /\ WF_vars(CompleteSharedBuild)
    /\ \A consumer \in Consumers: WF_vars(PublishTo(consumer))
    /\ WF_vars(RetireSharedResult)
    /\ WF_vars(CancelUnownedBuild)

ClosedConsumerHasNoPublication ==
    \A consumer \in Consumers:
        ~consumerOpen[consumer] => ~published[consumer]

CancellationRequiresNoConsumers ==
    producerState = "cancelled" => LiveConsumers = {}

SurvivingConsumerPreventsCancellation ==
    producerState = "building" /\ LiveConsumers # {} =>
        ~ENABLED CancelUnownedBuild

LiveDemandHasProgressWitness ==
    UnsatisfiedConsumers # {} =>
        \/ ENABLED SubmitSharedBuild
        \/ ENABLED CompleteSharedBuild
        \/ \E consumer \in Consumers: ENABLED PublishTo(consumer)

Ready ==
    /\ producerState \in {"idle", "cancelled"}
    /\ UnsatisfiedConsumers = {}
Constrained == FALSE
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyServedOrClosed ==
    [](UnsatisfiedConsumers # {} =>
       <>(UnsatisfiedConsumers = {} \/ LiveConsumers = {}))

=============================================================================
