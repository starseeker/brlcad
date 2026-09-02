---------------------- MODULE ObolResultAuthentication ----------------------
\* Authentication contract for asynchronous asset results.
\*
\* A dense occurrence index, stable route, or demand number is not sufficient
\* by itself.  Worker input and output carry both the semantic population epoch
\* and the demand revision which selected the work.  Replacement and demand
\* changes may overtake an in-flight result; the owner thread must retire that
\* result as superseded and submit a successor rather than applying it to the
\* new population.  This model is the prospective guard against the broader
\* family behind stale compact indices, allocator-address reuse, and late
\* callbacks crossing a scene rebuild.

EXTENDS Naturals, TLC

CONSTANTS MaxPopulationEpoch, MaxDemandRevision, MaxRouteRevision

RequestStates == {"idle", "inflight", "result"}
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
          acceptedPopulation,
          acceptedDemand,
          acceptedRoute,
          acceptedValid

vars == <<populationEpoch, demandRevision, routeRevision, requestState,
          requestPopulation, requestDemand, requestRoute, resultPopulation,
          resultDemand, resultRoute, acceptedPopulation, acceptedDemand,
          acceptedRoute, acceptedValid>>

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
    /\ acceptedPopulation \in 0..NoPopulation
    /\ acceptedDemand \in 0..NoDemand
    /\ acceptedRoute \in 0..NoRoute
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
    /\ acceptedPopulation = NoPopulation
    /\ acceptedDemand = NoDemand
    /\ acceptedRoute = NoRoute
    /\ acceptedValid = FALSE

SubmitCurrent ==
    /\ requestState = "idle"
    /\ ~acceptedValid
    /\ requestState' = "inflight"
    /\ requestPopulation' = populationEpoch
    /\ requestDemand' = demandRevision
    /\ requestRoute' = routeRevision
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    resultPopulation, resultDemand, resultRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute,
                    acceptedValid>>

CompleteResult ==
    /\ requestState = "inflight"
    /\ requestState' = "result"
    /\ resultPopulation' = requestPopulation
    /\ resultDemand' = requestDemand
    /\ resultRoute' = requestRoute
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute,
                    acceptedValid>>

ApplyCurrentResult ==
    /\ requestState = "result"
    /\ CurrentResult
    /\ requestState' = "idle"
    /\ acceptedPopulation' = resultPopulation
    /\ acceptedDemand' = resultDemand
    /\ acceptedRoute' = resultRoute
    /\ acceptedValid' = TRUE
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute>>

RejectSupersededResult ==
    /\ requestState = "result"
    /\ ~CurrentResult
    /\ requestState' = "idle"
    /\ UNCHANGED <<populationEpoch, demandRevision, routeRevision,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute,
                    acceptedValid>>

\* Population replacement is a semantic identity change even if storage
\* reuses the same dense slot or address.  The old presentation may remain as
\* a retained visual until a current commit, but it is no longer authenticated
\* as the current population.
ReplacePopulation ==
    /\ populationEpoch < MaxPopulationEpoch
    /\ populationEpoch' = populationEpoch + 1
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<demandRevision, routeRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute>>

ChangeDemand ==
    /\ demandRevision < MaxDemandRevision
    /\ demandRevision' = demandRevision + 1
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<populationEpoch, routeRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute>>

\* Source routing has an identity lifetime independent of both population and
\* requested quality.  A new provider or cache route may reuse the same dense
\* population slot and demand value while invalidating its predecessor's
\* result.
ChangeRoute ==
    /\ routeRevision < MaxRouteRevision
    /\ routeRevision' = routeRevision + 1
    /\ acceptedValid' = FALSE
    /\ UNCHANGED <<populationEpoch, demandRevision, requestState,
                    requestPopulation, requestDemand, requestRoute,
                    resultPopulation, resultDemand, resultRoute,
                    acceptedPopulation, acceptedDemand, acceptedRoute>>

Next ==
    \/ SubmitCurrent
    \/ CompleteResult
    \/ ApplyCurrentResult
    \/ RejectSupersededResult
    \/ ReplacePopulation
    \/ ChangeDemand
    \/ ChangeRoute

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(SubmitCurrent)
    /\ WF_vars(CompleteResult)
    /\ WF_vars(ApplyCurrentResult)
    /\ WF_vars(RejectSupersededResult)

AuthenticatedResultIsCurrent ==
    acceptedValid =>
        /\ acceptedPopulation = populationEpoch
        /\ acceptedDemand = demandRevision
        /\ acceptedRoute = routeRevision

StaleResultCannotApply ==
    requestState = "result" /\ ~CurrentResult =>
        /\ ~ENABLED ApplyCurrentResult
        /\ ENABLED RejectSupersededResult

CurrentNeedHasProgressWitness ==
    ~acceptedValid =>
        \/ ENABLED SubmitCurrent
        \/ ENABLED CompleteResult
        \/ ENABLED ApplyCurrentResult
        \/ ENABLED RejectSupersededResult

Ready == acceptedValid /\ requestState = "idle"
Constrained == FALSE
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyAuthenticatedAtFinalRevisions ==
    []((populationEpoch = MaxPopulationEpoch /\
        demandRevision = MaxDemandRevision /\
        routeRevision = MaxRouteRevision) => <>Ready)

=============================================================================
