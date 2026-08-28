--------------------------- MODULE ObolCapacitySearch --------------------------
\* Bounded renderer-capacity search below the progressive drawing pipeline.
\*
\* Candidate indices denote numeric allocation budgets for one immutable
\* revision tuple.  PopulationOf maps those budgets to complete, discrete,
\* visual-importance-ordered scene populations; adjacent budgets may select
\* the same PoP cuts.  Zero is the known-safe coverage population and
\* CandidateCount + 1 is an unsafe sentinel.  A previously unseen population
\* consumes at most SampleLimit measurements.  A different budget selecting
\* an already classified population reuses that result without repainting,
\* and either transition strictly narrows the safe/unsafe budget bracket.  A
\* quiet view searches first for the
\* preferred steady cadence.  If quality debt remains and the independent
\* static deadline is available, it may advance exactly once to the static
\* goal, preserving the steady-safe lower bound.  A ready pose-only view may
\* start at the static goal: its retained presentation is the initial
\* candidate and may be coarsened only if it misses that hard deadline.  No
\* timer, repaint, or throughput-EMA update may reopen either search.  A hard
\* deadline abort is an immediate unsafe classification of the active
\* candidate; it narrows this certificate instead of creating another one.
\* Candidate allocation and exact presentation are explicit producer
\* barriers.  Neither may consume timing samples, and measurement cannot
\* begin until both have completed.
\*
\* This model proves ownership and termination.  It deliberately does not
\* model the numeric frame classifier or perceptual ordering of candidates;
\* those require independent property tests and graphical qualification.

EXTENDS FiniteSets, Naturals, TLC

CONSTANT CandidateCount, SampleLimit, PopulationOf

ASSUME /\ CandidateCount > 0
       /\ SampleLimit > 0

Phases == {"choose", "allocate", "present", "measure", "goal_done", "terminal"}
Goals == {"steady", "static"}
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
          startAtStatic,
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
          goalTransitions

vars == <<phase, goal, trueSteadyCapacity, trueStaticCapacity, staticEligible,
          startAtStatic,
          safe, unsafe, candidate, samplesRemaining, measured,
          measuredSafePopulations, measuredUnsafePopulations,
          certificateRevision, priorWidth, narrowed, goalTransitions>>

Capacity(g) == IF g = "steady" THEN trueSteadyCapacity ELSE trueStaticCapacity
CandidateCapacity(g) ==
    Cardinality({budget \in Candidates : PopulationOf[budget] <= Capacity(g)})

TypeOK ==
    /\ phase \in Phases
    /\ goal \in Goals
    /\ trueSteadyCapacity \in 0..CandidateCount
    /\ trueStaticCapacity \in trueSteadyCapacity..CandidateCount
    /\ staticEligible \in BOOLEAN
    /\ startAtStatic \in BOOLEAN
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

BracketSound == safe <= CandidateCapacity(goal)
    /\ CandidateCapacity(goal) < unsafe

CandidateOwned ==
    phase \in {"allocate", "present", "measure"} =>
        /\ safe < candidate
        /\ candidate < unsafe
        /\ candidate \notin measured
        /\ IF PopulationOf[candidate] \in
                  (measuredSafePopulations \union measuredUnsafePopulations)
              THEN samplesRemaining = 0
              ELSE samplesRemaining > 0

TerminalCertificate ==
    phase = "terminal" =>
        /\ unsafe = safe + 1
        /\ safe = CandidateCapacity(goal)
        /\ certificateRevision = 1
        /\ samplesRemaining = 0
        /\ \/ goal = "static"
           \/ ~staticEligible
           \/ safe = CandidateCount

GoalMonotonic ==
    /\ (goal = "steady" => /\ ~startAtStatic
                              /\ goalTransitions = 0)
    /\ (goal = "static" => goalTransitions = 1)

StrictNarrowing == narrowed => unsafe - safe < priorWidth

Init ==
    /\ phase = "choose"
    /\ startAtStatic \in BOOLEAN
    /\ goal = IF startAtStatic THEN "static" ELSE "steady"
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
    /\ goalTransitions = IF startAtStatic THEN 1 ELSE 0

ChooseCandidate ==
    /\ phase = "choose"
    /\ unsafe > safe + 1
    /\ candidate' = safe + (unsafe - safe) \div 2
    /\ samplesRemaining' =
          IF PopulationOf[safe + (unsafe - safe) \div 2] \in
                  (measuredSafePopulations \union measuredUnsafePopulations)
              THEN 0 ELSE SampleLimit
    /\ phase' = "allocate"
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    certificateRevision, priorWidth, goalTransitions>>

ApplyCandidateAllocation ==
    /\ phase = "allocate"
    /\ phase' = "present"
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions>>

PresentExactCandidate ==
    /\ phase = "present"
    /\ phase' = "measure"
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, certificateRevision,
                    priorWidth, goalTransitions>>

ConsumeSample ==
    /\ phase = "measure"
    /\ samplesRemaining > 1
    /\ samplesRemaining' = samplesRemaining - 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<phase, goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, candidate, measured,
                    measuredSafePopulations, measuredUnsafePopulations,
                    certificateRevision, priorWidth, goalTransitions>>

ClassifyCandidate ==
    /\ phase = "measure"
    /\ samplesRemaining = 1
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
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ narrowed' = TRUE
    /\ IF unsafe' = safe' + 1
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, goalTransitions>>

RejectAtDeadline ==
    /\ phase = "measure"
    /\ samplesRemaining > 0
    /\ PopulationOf[candidate] > Capacity(goal)
    /\ priorWidth' = unsafe - safe
    /\ safe' = safe
    /\ unsafe' = candidate
    /\ measured' = measured \cup {candidate}
    /\ measuredSafePopulations' = measuredSafePopulations
    /\ measuredUnsafePopulations' =
          measuredUnsafePopulations \union {PopulationOf[candidate]}
    /\ candidate' = 0
    /\ samplesRemaining' = 0
    /\ narrowed' = TRUE
    /\ IF unsafe' = safe' + 1
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, goalTransitions>>

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
    /\ candidate' = 0
    /\ narrowed' = TRUE
    /\ IF unsafe' = safe' + 1
          THEN /\ phase' = "goal_done"
               /\ certificateRevision' = 0
          ELSE /\ phase' = "choose"
               /\ certificateRevision' = 0
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, samplesRemaining,
                    measuredSafePopulations, measuredUnsafePopulations,
                    goalTransitions>>

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
    /\ UNCHANGED <<trueSteadyCapacity, trueStaticCapacity, staticEligible,
                    startAtStatic,
                    safe>>

PublishTerminalCertificate ==
    /\ phase = "goal_done"
    /\ ~(/\ goal = "steady"
          /\ staticEligible
          /\ safe < CandidateCount)
    /\ phase' = "terminal"
    /\ certificateRevision' = 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, candidate,
                    samplesRemaining, measured, measuredSafePopulations,
                    measuredUnsafePopulations, priorWidth, goalTransitions>>

Next ==
    \/ ChooseCandidate
    \/ ApplyCandidateAllocation
    \/ PresentExactCandidate
    \/ ConsumeSample
    \/ ClassifyCandidate
    \/ RejectAtDeadline
    \/ ReusePopulation
    \/ AdvanceToStaticGoal
    \/ PublishTerminalCertificate

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ChooseCandidate)
    /\ WF_vars(ApplyCandidateAllocation)
    /\ WF_vars(PresentExactCandidate)
    /\ WF_vars(ConsumeSample)
    /\ WF_vars(ClassifyCandidate)
    /\ WF_vars(ReusePopulation)
    /\ WF_vars(AdvanceToStaticGoal)
    /\ WF_vars(PublishTerminalCertificate)

EventuallyCertified == <> (phase = "terminal")

=============================================================================
