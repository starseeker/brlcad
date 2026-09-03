---------------------- MODULE ObolResultAuthentication ----------------------
\* Authentication and typed-disposition contract for asynchronous results.
\*
\* A dense occurrence index, stable route, or demand number is not sufficient
\* by itself.  Worker input and output carry both the semantic population epoch
\* and the demand revision which selected the work.  Replacement and demand
\* changes may overtake an in-flight result; the owner thread must retire that
\* result as superseded and submit a successor rather than applying it to the
\* new population or recording its terminal error against the successor.
\*
\* Provider statuses are abstracted to the three owner-thread dispositions:
\* publish useful output, record an exact-current terminal failure, or retry
\* current demand.  MaxRetryResults bounds retryable completions only to keep
\* the liveness state space finite; it does not bound production retries.

EXTENDS Naturals, TLC

CONSTANTS MaxPopulationEpoch, MaxDemandRevision, MaxRouteRevision,
          MaxRetryResults

RequestStates == {"idle", "inflight", "result"}
ResultOutcomes == {"publish", "terminalFailure", "retry"}
AcceptedOutcomes == {"none", "published", "failed"}
NoPopulation == MaxPopulationEpoch + 1
NoDemand == MaxDemandRevision + 1
NoRoute == MaxRouteRevision + 1

VARIABLES populationEpoch,
          demandRevision,
          routeRevision,
          requestState,
          requestPopulation,
          requestDemand,
          requestRoute,
          resultPopulation,
          resultDemand,
          resultRoute,
          resultOutcome,
          retryResultsRemaining,
          acceptedPopulation,
          acceptedDemand,
          acceptedRoute,
          acceptedOutcome,
          acceptedValid

vars == <<populationEpoch, demandRevision, routeRevision, requestState,
          requestPopulation, requestDemand, requestRoute, resultPopulation,
          resultDemand, resultRoute, resultOutcome, retryResultsRemaining,
          acceptedPopulation, acceptedDemand, acceptedRoute, acceptedOutcome,
          acceptedValid>>

CurrentResult ==
    /\ resultPopulation = populationEpoch
    /\ resultDemand = demandRevision
    /\ resultRoute = routeRevision

TypeOK ==
    /\ populationEpoch \in 0..MaxPopulationEpoch
    /\ demandRevision \in 0..MaxDemandRevision
    /\ routeRevision \in 0..MaxRouteRevision
    /\ requestState \in RequestStates
    /\ requestPopulation \in 0..MaxPopulationEpoch
    /\ requestDemand \in 0..MaxDemandRevision
    /\ requestRoute \in 0..MaxRouteRevision
    /\ resultPopulation \in 0..NoPopulation
    /\ resultDemand \in 0..NoDemand
    /\ resultRoute \in 0..NoRoute
    /\ resultOutcome \in ResultOutcomes \cup {"none"}
    /\ retryResultsRemaining \in 0..MaxRetryResults
    /\ acceptedPopulation \in 0..NoPopulation
    /\ acceptedDemand \in 0..NoDemand
    /\ acceptedRoute \in 0..NoRoute
    /\ acceptedOutcome \in AcceptedOutcomes
    /\ acceptedValid \in BOOLEAN

Init ==
    /\ populationEpoch = 0
    /\ demandRevision = 0
    /\ routeRevision = 0
    /\ requestState = "idle"
    /\ requestPopulation = 0
    /\ requestDemand = 0
    /\ requestRoute = 0
    /\ resultPopulation = NoPopulation
    /\ resultDemand = NoDemand
    /\ resultRoute = NoRoute
    /\ resultOutcome = "none"
    /\ retryResultsRemaining = MaxRetryResults
    /\ acceptedPopulation = NoPopulation
    /\ acceptedDemand = NoDemand
    /\ acceptedRoute = NoRoute
    /\ acceptedOutcome = "none"
    /\ acceptedValid = FALSE

SubmitCurrent ==
    /\ requestState = "idle"
    /\ ~acceptedValid
    /\ requestState' = "inflight"
    /\ requestPopulation' = populationEpoch
    /\ requestDemand' = demandRevision
    /\ requestRoute' = routeRevision
    /\ resultOutcome' = "none"
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    resultPopulation, resultDemand, resultRoute,
                    retryResultsRemaining, acceptedPopulation,
                    acceptedDemand, acceptedRoute, acceptedOutcome,
                    acceptedValid>>

CompleteResult(outcome) ==
    /\ requestState = "inflight"
    /\ outcome \in ResultOutcomes
    /\ outcome = "retry" => retryResultsRemaining > 0
    /\ requestState' = "result"
    /\ resultPopulation' = requestPopulation
    /\ resultDemand' = requestDemand
    /\ resultRoute' = requestRoute
    /\ resultOutcome' = outcome
    /\ retryResultsRemaining' =
        IF outcome = "retry"
        THEN retryResultsRemaining - 1
        ELSE retryResultsRemaining
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute,
                    acceptedOutcome, acceptedValid>>

CompleteAnyResult ==
    \E outcome \in ResultOutcomes: CompleteResult(outcome)

ApplyCurrentResult ==
    /\ requestState = "result"
    /\ CurrentResult
    /\ resultOutcome = "publish"
    /\ requestState' = "idle"
    /\ resultOutcome' = "none"
    /\ acceptedPopulation' = resultPopulation
    /\ acceptedDemand' = resultDemand
    /\ acceptedRoute' = resultRoute
    /\ acceptedOutcome' = "published"
    /\ acceptedValid' = TRUE
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    retryResultsRemaining>>

RecordCurrentFailure ==
    /\ requestState = "result"
    /\ CurrentResult
    /\ resultOutcome = "terminalFailure"
    /\ requestState' = "idle"
    /\ resultOutcome' = "none"
    /\ acceptedPopulation' = resultPopulation
    /\ acceptedDemand' = resultDemand
    /\ acceptedRoute' = resultRoute
    /\ acceptedOutcome' = "failed"
    /\ acceptedValid' = TRUE
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    retryResultsRemaining>>

RetryCurrentResult ==
    /\ requestState = "result"
    /\ CurrentResult
    /\ resultOutcome = "retry"
    /\ requestState' = "idle"
    /\ resultOutcome' = "none"
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    retryResultsRemaining, acceptedPopulation,
                    acceptedDemand, acceptedRoute, acceptedOutcome,
                    acceptedValid>>

RejectSupersededResult ==
    /\ requestState = "result"
    /\ ~CurrentResult
    /\ requestState' = "idle"
    /\ resultOutcome' = "none"
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    retryResultsRemaining, acceptedPopulation,
                    acceptedDemand, acceptedRoute, acceptedOutcome,
                    acceptedValid>>

\* Population replacement is a semantic identity change even if storage
\* reuses the same dense slot or address.  The old presentation may remain as
\* a retained visual until a current commit, but it is no longer authenticated
\* as the current population.
ReplacePopulation ==
    /\ populationEpoch < MaxPopulationEpoch
    /\ populationEpoch' = populationEpoch + 1
    /\ acceptedOutcome' = "none"
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<demandRevision, routeRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute, resultOutcome,
                    retryResultsRemaining, acceptedPopulation, acceptedDemand,
                    acceptedRoute>>

ChangeDemand ==
    /\ demandRevision < MaxDemandRevision
    /\ demandRevision' = demandRevision + 1
    /\ acceptedOutcome' = "none"
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<populationEpoch, routeRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute, resultOutcome,
                    retryResultsRemaining, acceptedPopulation, acceptedDemand,
                    acceptedRoute>>

\* Source routing has an identity lifetime independent of both population and
\* requested quality.  A new provider or cache route may reuse the same dense
\* population slot and demand value while invalidating its predecessor's
\* result.
ChangeRoute ==
    /\ routeRevision < MaxRouteRevision
    /\ routeRevision' = routeRevision + 1
    /\ acceptedOutcome' = "none"
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<populationEpoch, demandRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute, resultOutcome,
                    retryResultsRemaining, acceptedPopulation, acceptedDemand,
                    acceptedRoute>>

Next ==
    \/ SubmitCurrent
    \/ CompleteAnyResult
    \/ ApplyCurrentResult
    \/ RecordCurrentFailure
    \/ RetryCurrentResult
    \/ RejectSupersededResult
    \/ ReplacePopulation
    \/ ChangeDemand
    \/ ChangeRoute

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SubmitCurrent)
    /\ WF_vars(CompleteAnyResult)
    /\ WF_vars(ApplyCurrentResult)
    /\ WF_vars(RecordCurrentFailure)
    /\ WF_vars(RetryCurrentResult)
    /\ WF_vars(RejectSupersededResult)

AuthenticatedResultIsCurrent ==
    acceptedValid =>
        /\ acceptedPopulation = populationEpoch
        /\ acceptedDemand = demandRevision
        /\ acceptedRoute = routeRevision

AcceptedOutcomeIsTyped ==
    acceptedValid = (acceptedOutcome \in {"published", "failed"})

StaleResultCannotApply ==
    requestState = "result" /\ ~CurrentResult =>
        /\ ~ENABLED ApplyCurrentResult
        /\ ~ENABLED RecordCurrentFailure
        /\ ENABLED RejectSupersededResult

SupersededResultCreatesNoOutcome ==
    requestState = "result" /\ ~CurrentResult => ~acceptedValid

CurrentNeedHasProgressWitness ==
    ~acceptedValid =>
        \/ ENABLED SubmitCurrent
        \/ ENABLED CompleteAnyResult
        \/ ENABLED ApplyCurrentResult
        \/ ENABLED RecordCurrentFailure
        \/ ENABLED RetryCurrentResult
        \/ ENABLED RejectSupersededResult

Ready ==
    /\ acceptedValid
    /\ acceptedOutcome = "published"
    /\ requestState = "idle"
Constrained == FALSE
Failed ==
    /\ acceptedValid
    /\ acceptedOutcome = "failed"
    /\ requestState = "idle"
Terminal == Ready \/ Constrained \/ Failed

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyAuthenticatedAtFinalRevisions ==
    []((populationEpoch = MaxPopulationEpoch /\
        demandRevision = MaxDemandRevision /\
        routeRevision = MaxRouteRevision) => <>Terminal)

=============================================================================
