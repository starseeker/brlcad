--------------------------- MODULE ObolCapacitySearch --------------------------
\* Bounded renderer-capacity search below the progressive drawing pipeline.
\*
\* Candidate indices denote allocator-realizable budgets for one immutable
\* revision tuple.  That tuple includes the aggregate point/box
\* classification rule and pixel threshold; changing either starts a new
\* instance of this machine.  The candidate domain spans the complete
\* currently resident pixel-demand cost.  PopulationOf maps those valid scene
\* budgets to complete, discrete,
\* visual-importance-ordered scene populations; adjacent budgets may select
\* the same PoP cuts.  Zero is the known-safe coverage population and
\* CandidateCount + 1 is an unsafe sentinel.  A previously unseen population
\* consumes at most SampleLimit measurements.  A different budget selecting
\* an already classified population reuses that result without repainting,
\* and either transition strictly narrows the safe/unsafe budget bracket.  A
\* quiet view searches first for the
\* preferred steady cadence.  If quality debt remains and the independent
\* static deadline is available, it may advance exactly once to the static
\* goal, preserving the steady-safe lower bound.  Any ready view with a
\* retained mesh presentation, after pose or zoom, starts at the static goal:
\* its retained presentation is the initial candidate and may be coarsened
\* only if it misses that hard deadline.  Zoom still recomputes pixel demand;
\* this rule governs the capacity lower bound, not occurrence-cut reuse.  No
\* timer, repaint, or throughput-EMA update may reopen either search.  A hard
\* deadline abort is an immediate unsafe classification of exactly the active
\* candidate; it narrows this certificate instead of creating another one.
\* A throughput-derived recovery hint may select a lower successor, but may
\* not move the unsafe boundary below the candidate which actually missed.
\* That distinction is essential when adjacent budgets map to different,
\* widely separated PoP populations.
\* SearchMinimum is the allocator's irreducible lower endpoint.  When a
\* protected visual floor is active, production sets it to that floor's exact
\* cost rather than allowing the allocator to select the floor as an override
\* for a smaller capacity candidate.  Consequently every selected population
\* remains inside its candidate budget and a terminal certificate cannot be
\* followed by restoration of an out-of-domain floor.
\* Capacity is a responsiveness guard, not a visual-quality objective.  Once
\* the final goal has both a proven-safe population and a rejected richer
\* population, it may publish one midpoint bridge before settling.  This
\* recovers a wide unused bracket without exposing an open-ended binary search
\* or repeated quality staircase.  CandidateAttemptLimit
\* independently bounds cold recovery before any safe population is known.
\* Candidate allocation and exact presentation are explicit producer
\* barriers.  Neither may consume timing samples, and measurement cannot
\* begin until both have completed.  Inexact presentations and unusable
\* measurement frames share one bounded consecutive-attempt certificate; a
\* valid frame resets it and no retry is an unranked stutter.  An older
\* population-presentation barrier
\* may still be pending when a completed sample chooses the next allocation.
\* That barrier owns the already completed frame first, then transfers the
\* exact candidate allocation without retaining presentation debt.
\*
\* This model proves ownership, bounded candidate publication, one-sided
\* convergence after a safe population exists, and termination.  Its monotone
\* PopulationOf
\* assumption is discharged by ObolRetainedAllocationPrefix.tla and the
\* exhaustive retained-allocation C++ budget sweep.  It deliberately does not
\* model the numeric frame classifier or perceptual ordering arithmetic;
\* those require independent property tests and graphical qualification.

EXTENDS FiniteSets, Naturals, TLC

CONSTANT CandidateCount, SearchMinimum, SampleLimit,
         CandidateAttemptLimit, InvalidFrameAttemptLimit, PopulationOf

ASSUME /\ CandidateCount > 0
       /\ SearchMinimum \in 1..CandidateCount
       /\ SampleLimit > 0
       /\ CandidateAttemptLimit > 0
       /\ InvalidFrameAttemptLimit > 0

Phases == {"choose", "allocate", "present", "measure", "goal_done",
           "terminal", "unmeasurable"}
Goals == {"steady", "static"}
HandoffKinds == {"cold", "readyPose", "readyZoom"}
Candidates == 1..CandidateCount
Populations == 1..CandidateCount
PairedPopulationMap ==
    [candidate \in Candidates |-> (candidate + 1) \div 2]

ASSUME /\ PopulationOf \in [Candidates -> Populations]
       /\ \A candidate \in Candidates : PopulationOf[candidate] <= candidate
       /\ \A lower, upper \in Candidates :
              lower <= upper => PopulationOf[lower] <= PopulationOf[upper]

VARIABLES phase,
          goal,
          trueSteadyCapacity,
          trueStaticCapacity,
          staticEligible,
          handoffKind,
          safe,
          unsafe,
          candidate,
          samplesRemaining,
          measured,
          measuredSafePopulations,
          measuredUnsafePopulations,
          certificateRevision,
          priorWidth,
          narrowed,
          goalTransitions,
          goalAttempts,
          candidatePublications,
          postSafeFailure,
          bracketBridgeUsed,
          priorFrameBarrier,
          revisionHasPlan,
          invalidFrameAttempts

vars == <<phase, goal, trueSteadyCapacity, trueStaticCapacity, staticEligible,
          handoffKind,
          safe, unsafe, candidate, samplesRemaining, measured,
          measuredSafePopulations, measuredUnsafePopulations,
          certificateRevision, priorWidth, narrowed, goalTransitions,
          goalAttempts, candidatePublications, postSafeFailure,
          bracketBridgeUsed, priorFrameBarrier, revisionHasPlan,
          invalidFrameAttempts>>

Capacity(g) == IF g = "steady" THEN trueSteadyCapacity ELSE trueStaticCapacity
StartAtStatic == handoffKind \in {"readyPose", "readyZoom"}
CandidateCapacity(g) ==
    Cardinality({budget \in Candidates : PopulationOf[budget] <= Capacity(g)})

TypeOK ==
    /\ phase \in Phases
    /\ goal \in Goals
    /\ trueSteadyCapacity \in 0..CandidateCount
    /\ trueStaticCapacity \in trueSteadyCapacity..CandidateCount
    /\ staticEligible \in BOOLEAN
    /\ handoffKind \in HandoffKinds
    /\ safe \in 0..CandidateCount
    /\ unsafe \in 1..(CandidateCount + 1)
    /\ candidate \in 0..CandidateCount
    /\ samplesRemaining \in 0..SampleLimit
    /\ measured \subseteq Candidates
    /\ measuredSafePopulations \subseteq Populations
    /\ measuredUnsafePopulations \subseteq Populations
    /\ measuredSafePopulations \cap measuredUnsafePopulations = {}
    /\ certificateRevision \in 0..1
    /\ priorWidth \in 1..(CandidateCount + 1)
    /\ narrowed \in BOOLEAN
    /\ goalTransitions \in 0..1
    /\ goalAttempts \in 0..CandidateAttemptLimit
    /\ candidatePublications \in 0..(2 * CandidateAttemptLimit)
    /\ postSafeFailure \in BOOLEAN
    /\ bracketBridgeUsed \in BOOLEAN
    /\ priorFrameBarrier \in BOOLEAN
    /\ revisionHasPlan \in BOOLEAN
    /\ invalidFrameAttempts \in 0..InvalidFrameAttemptLimit

BracketSound == safe <= CandidateCapacity(goal)
    /\ CandidateCapacity(goal) < unsafe

CandidateOwned ==
    phase \in {"allocate", "present", "measure"} =>
        /\ safe < candidate
        /\ candidate < unsafe
        /\ candidate >= SearchMinimum
        /\ candidate \notin measured
        /\ IF PopulationOf[candidate] \in
                  (measuredSafePopulations \union measuredUnsafePopulations)
              THEN samplesRemaining = 0
              ELSE samplesRemaining > 0

PriorFrameBarrierOwnsAllocation ==
    priorFrameBarrier => phase = "allocate"

PlanRevisionOwnership ==
    /\ (phase = "allocate" => revisionHasPlan = priorFrameBarrier)
    /\ (phase \in {"present", "measure"} => revisionHasPlan)

TerminalCertificate ==
    phase = "terminal" =>
        /\ certificateRevision = 1
        /\ samplesRemaining = 0
        /\ safe <= CandidateCapacity(goal)
        /\ CandidateCapacity(goal) < unsafe
        /\ \/ unsafe = safe + 1
           \/ safe = CandidateCount
           \/ goalAttempts = CandidateAttemptLimit
           \/ /\ safe < SearchMinimum
              /\ unsafe <= SearchMinimum
           \/ /\ safe > 0
              /\ unsafe <= CandidateCount
        /\ \/ goal = "static"
           \/ ~staticEligible
           \/ safe = CandidateCount

UnmeasurableCertificate ==
    phase = "unmeasurable" =>
        /\ certificateRevision = 1
        /\ samplesRemaining = 0

GoalMonotonic ==
    /\ (goal = "steady" => /\ ~StartAtStatic
                              /\ goalTransitions = 0)
    /\ (goal = "static" => goalTransitions = 1)

RetainedHandoffStartsAtStatic == StartAtStatic => goal = "static"

CandidatePublicationBound ==
    candidatePublications <=
        IF StartAtStatic THEN CandidateAttemptLimit
        ELSE 2 * CandidateAttemptLimit

PostSafeFailureTransfersOnce ==
    postSafeFailure =>
        \/ phase \in {"goal_done", "terminal"}
        \/ /\ phase = "choose"
           /\ ~bracketBridgeUsed

StrictNarrowing == narrowed => unsafe - safe < priorWidth

Init ==
    /\ phase = "choose"
    /\ handoffKind \in HandoffKinds
    /\ goal = IF StartAtStatic THEN "static" ELSE "steady"
    /\ trueSteadyCapacity \in 0..CandidateCount
    /\ trueStaticCapacity \in trueSteadyCapacity..CandidateCount
    /\ staticEligible \in BOOLEAN
    /\ safe = 0
    /\ unsafe = CandidateCount + 1
    /\ candidate = 0
    /\ samplesRemaining = 0
    /\ measured = {}
    /\ measuredSafePopulations = {}
    /\ measuredUnsafePopulations = {}
    /\ certificateRevision = 0
    /\ priorWidth = CandidateCount + 1
    /\ narrowed = FALSE
    /\ goalTransitions = IF StartAtStatic THEN 1 ELSE 0
    /\ goalAttempts = 0
    /\ candidatePublications = 0
    /\ postSafeFailure = FALSE
    /\ bracketBridgeUsed = FALSE
    /\ priorFrameBarrier = FALSE
    /\ revisionHasPlan = FALSE
    /\ invalidFrameAttempts = 0

ChooseCandidate ==
    /\ phase = "choose"
    /\ unsafe > IF safe + 1 > SearchMinimum
                  THEN safe + 1 ELSE SearchMinimum
    /\ goalAttempts < CandidateAttemptLimit
    /\ candidate' =
          IF safe + (unsafe - safe) \div 2 < SearchMinimum
          THEN SearchMinimum
          ELSE safe + (unsafe - safe) \div 2
    /\ samplesRemaining' =
          IF PopulationOf[candidate'] \in
                  (measuredSafePopulations \union measuredUnsafePopulations)
              THEN 0 ELSE SampleLimit
    /\ \E barrier \in BOOLEAN :
          /\ priorFrameBarrier' = barrier
          /\ revisionHasPlan' = barrier
    /\ phase' = "allocate"
    /\ goalAttempts' = goalAttempts + 1
    /\ candidatePublications' = candidatePublications + 1
    /\ postSafeFailure' = FALSE
    /\ bracketBridgeUsed' =
          (bracketBridgeUsed \/ (safe > 0 /\ unsafe <= CandidateCount))
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    certificateRevision, priorWidth, goalTransitions,
                    invalidFrameAttempts>>

ConsumePriorFrameBarrier ==
    /\ phase = "allocate"
    /\ priorFrameBarrier
    /\ priorFrameBarrier' = FALSE
    /\ revisionHasPlan' = FALSE
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, phase,
                    bracketBridgeUsed, invalidFrameAttempts>>

ApplyCandidateAllocation ==
    /\ phase = "allocate"
    /\ ~priorFrameBarrier
    /\ phase' = "present"
    /\ invalidFrameAttempts' = 0
    /\ revisionHasPlan' = TRUE
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure,
                    bracketBridgeUsed, priorFrameBarrier>>

PresentExactCandidate ==
    /\ phase = "present"
    /\ phase' = "measure"
    /\ invalidFrameAttempts' = 0
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, priorFrameBarrier,
                    revisionHasPlan, bracketBridgeUsed>>

RecordInexactPresentation ==
    /\ phase = "present"
    /\ invalidFrameAttempts < InvalidFrameAttemptLimit
    /\ invalidFrameAttempts' = invalidFrameAttempts + 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<phase, goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, priorFrameBarrier,
                    revisionHasPlan, bracketBridgeUsed>>

RetireUnmeasurablePresentation ==
    /\ phase = "present"
    /\ invalidFrameAttempts = InvalidFrameAttemptLimit
    /\ phase' = "unmeasurable"
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ certificateRevision' = 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, priorFrameBarrier,
                    revisionHasPlan, invalidFrameAttempts,
                    bracketBridgeUsed>>

ResolvePresentation ==
    PresentExactCandidate \/ RecordInexactPresentation
        \/ RetireUnmeasurablePresentation

RecordInvalidMeasurement ==
    /\ phase = "measure"
    /\ invalidFrameAttempts < InvalidFrameAttemptLimit
    /\ invalidFrameAttempts' = invalidFrameAttempts + 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<phase, goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, priorFrameBarrier,
                    revisionHasPlan, bracketBridgeUsed>>

RetireUnmeasurableMeasurement ==
    /\ phase = "measure"
    /\ invalidFrameAttempts = InvalidFrameAttemptLimit
    /\ phase' = "unmeasurable"
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ certificateRevision' = 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    priorWidth, goalTransitions, goalAttempts,
                    candidatePublications, postSafeFailure, priorFrameBarrier,
                    revisionHasPlan, invalidFrameAttempts,
                    bracketBridgeUsed>>

ResolveInvalidMeasurement ==
    RecordInvalidMeasurement \/ RetireUnmeasurableMeasurement

ConsumeSample ==
    /\ phase = "measure"
    /\ samplesRemaining > 1
    /\ samplesRemaining' = samplesRemaining - 1
    /\ invalidFrameAttempts' = 0
    /\ narrowed' = FALSE
    /\ UNCHANGED <<phase, goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    certificateRevision, priorWidth, goalTransitions,
                    goalAttempts, candidatePublications, postSafeFailure,
                    priorFrameBarrier, revisionHasPlan, bracketBridgeUsed>>

ClassifyCandidate ==
    /\ phase = "measure"
    /\ samplesRemaining = 1
    /\ ~( /\ goal = "steady"
           /\ staticEligible
           /\ PopulationOf[candidate] > trueSteadyCapacity
           /\ PopulationOf[candidate] <= trueStaticCapacity )
    /\ priorWidth' = unsafe - safe
    /\ measured' = measured \cup {candidate}
    /\ IF PopulationOf[candidate] <= Capacity(goal)
          THEN /\ safe' = candidate
               /\ unsafe' = unsafe
               /\ measuredSafePopulations' =
                    measuredSafePopulations \union {PopulationOf[candidate]}
               /\ measuredUnsafePopulations' = measuredUnsafePopulations
          ELSE /\ safe' = safe
               /\ unsafe' = candidate
               /\ measuredSafePopulations' = measuredSafePopulations
               /\ measuredUnsafePopulations' =
                    measuredUnsafePopulations \union {PopulationOf[candidate]}
    /\ postSafeFailure' =
          (safe > 0 /\ PopulationOf[candidate] > Capacity(goal))
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ invalidFrameAttempts' = 0
    /\ narrowed' = TRUE
    /\ IF \/ unsafe' = safe' + 1
           \/ safe' = CandidateCount
           \/ goalAttempts = CandidateAttemptLimit
           \/ /\ safe' < SearchMinimum
              /\ unsafe' <= SearchMinimum
           \/ /\ safe' > 0
              /\ unsafe' <= CandidateCount
              /\ \/ bracketBridgeUsed
                 \/ /\ goal = "steady"
                    /\ staticEligible
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, goalTransitions,
                    goalAttempts, candidatePublications,
                    priorFrameBarrier, revisionHasPlan, bracketBridgeUsed>>

PromoteStaticSafeCandidate ==
    /\ phase = "measure"
    /\ samplesRemaining = 1
    /\ goal = "steady"
    /\ staticEligible
    /\ PopulationOf[candidate] > trueSteadyCapacity
    /\ PopulationOf[candidate] <= trueStaticCapacity
    /\ goal' = "static"
    /\ safe' = candidate
    /\ unsafe' = CandidateCount + 1
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ measured' = {}
    /\ measuredSafePopulations' = {}
    /\ measuredUnsafePopulations' = {}
    /\ certificateRevision' = 0
    /\ priorWidth' = CandidateCount + 1 - candidate
    /\ narrowed' = FALSE
    /\ goalTransitions' = 1
    /\ goalAttempts' = 0
    /\ postSafeFailure' = FALSE
    /\ bracketBridgeUsed' = FALSE
    /\ invalidFrameAttempts' = 0
    /\ phase' = IF candidate = CandidateCount THEN "goal_done" ELSE "choose"
    /\ UNCHANGED <<trueSteadyCapacity, trueStaticCapacity, handoffKind,
                    staticEligible, candidatePublications,
                    priorFrameBarrier, revisionHasPlan>>

RejectAtDeadline ==
    /\ phase \in {"allocate", "present", "measure"}
    /\ ~priorFrameBarrier
    /\ PopulationOf[candidate] > Capacity(goal)
    /\ priorWidth' = unsafe - safe
    /\ safe' = safe
    /\ unsafe' = candidate
    /\ measured' = measured \cup {candidate}
    /\ measuredSafePopulations' = measuredSafePopulations
    /\ measuredUnsafePopulations' =
          measuredUnsafePopulations \union {PopulationOf[candidate]}
    /\ postSafeFailure' = (safe > 0)
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ narrowed' = TRUE
    /\ IF \/ unsafe' = safe' + 1
           \/ goalAttempts = CandidateAttemptLimit
           \/ /\ safe' < SearchMinimum
              /\ unsafe' <= SearchMinimum
           \/ /\ safe' > 0
              /\ unsafe' <= CandidateCount
              /\ \/ bracketBridgeUsed
                 \/ /\ goal = "steady"
                    /\ staticEligible
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, goalTransitions,
                    goalAttempts, candidatePublications,
                    priorFrameBarrier, revisionHasPlan,
                    invalidFrameAttempts, bracketBridgeUsed>>

ReusePopulation ==
    /\ phase = "measure"
    /\ samplesRemaining = 0
    /\ PopulationOf[candidate] \in
          (measuredSafePopulations \union measuredUnsafePopulations)
    /\ priorWidth' = unsafe - safe
    /\ measured' = measured \union {candidate}
    /\ IF PopulationOf[candidate] \in measuredSafePopulations
          THEN /\ safe' = candidate
               /\ unsafe' = unsafe
          ELSE /\ safe' = safe
               /\ unsafe' = candidate
    /\ postSafeFailure' =
          (safe > 0 /\ PopulationOf[candidate] \in
              measuredUnsafePopulations)
    /\ candidate' = 0
    /\ narrowed' = TRUE
    /\ IF \/ unsafe' = safe' + 1
           \/ safe' = CandidateCount
           \/ goalAttempts = CandidateAttemptLimit
           \/ /\ safe' < SearchMinimum
              /\ unsafe' <= SearchMinimum
           \/ /\ safe' > 0
              /\ unsafe' <= CandidateCount
              /\ \/ bracketBridgeUsed
                 \/ /\ goal = "steady"
                    /\ staticEligible
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, samplesRemaining,
                    measuredSafePopulations, measuredUnsafePopulations,
                    goalTransitions, goalAttempts, candidatePublications,
                    priorFrameBarrier, revisionHasPlan,
                    invalidFrameAttempts, bracketBridgeUsed>>

AdvanceToStaticGoal ==
    /\ phase = "goal_done"
    /\ goal = "steady"
    /\ staticEligible
    /\ safe < CandidateCount
    /\ goal' = "static"
    /\ phase' = "choose"
    /\ unsafe' = CandidateCount + 1
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ measured' = {}
    /\ measuredSafePopulations' = {}
    /\ measuredUnsafePopulations' = {}
    /\ certificateRevision' = 0
    /\ priorWidth' = CandidateCount + 1 - safe
    /\ narrowed' = FALSE
    /\ goalTransitions' = 1
    /\ goalAttempts' = 0
    /\ postSafeFailure' = FALSE
    /\ bracketBridgeUsed' = FALSE
    /\ UNCHANGED <<trueSteadyCapacity, trueStaticCapacity, staticEligible,
                    handoffKind,
                    safe, candidatePublications,
                    priorFrameBarrier, revisionHasPlan,
                    invalidFrameAttempts>>

PublishTerminalCertificate ==
    /\ phase = "goal_done"
    /\ ~(/\ goal = "steady"
          /\ staticEligible
          /\ safe < CandidateCount)
    /\ phase' = "terminal"
    /\ certificateRevision' = 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, handoffKind, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, priorWidth, goalTransitions,
                    goalAttempts, candidatePublications, postSafeFailure,
                    priorFrameBarrier, revisionHasPlan,
                    invalidFrameAttempts, bracketBridgeUsed>>

Next ==
    \/ ChooseCandidate
    \/ ConsumePriorFrameBarrier
    \/ ApplyCandidateAllocation
    \/ ResolvePresentation
    \/ ResolveInvalidMeasurement
    \/ ConsumeSample
    \/ PromoteStaticSafeCandidate
    \/ ClassifyCandidate
    \/ RejectAtDeadline
    \/ ReusePopulation
    \/ AdvanceToStaticGoal
    \/ PublishTerminalCertificate

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ChooseCandidate)
    /\ WF_vars(ConsumePriorFrameBarrier)
    /\ WF_vars(ApplyCandidateAllocation)
    /\ WF_vars(ResolvePresentation)
    /\ WF_vars(ResolveInvalidMeasurement)
    /\ WF_vars(ConsumeSample)
    /\ WF_vars(PromoteStaticSafeCandidate)
    /\ WF_vars(ClassifyCandidate)
    /\ WF_vars(ReusePopulation)
    /\ WF_vars(AdvanceToStaticGoal)
    /\ WF_vars(PublishTerminalCertificate)

EventuallySettled == <> (phase \in {"terminal", "unmeasurable"})

DeadlockOnlyAtTerminal ==
    ~ENABLED <<Next>>_vars => phase \in {"terminal", "unmeasurable"}

=============================================================================
