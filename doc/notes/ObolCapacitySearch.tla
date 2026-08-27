--------------------------- MODULE ObolCapacitySearch --------------------------
\* Bounded renderer-capacity search below the progressive drawing pipeline.
\*
\* Candidate indices denote complete, visual-importance-ordered scene
\* allocations for one immutable revision tuple.  Zero is the known-safe
\* coverage population and CandidateCount + 1 is an unsafe sentinel.  A
\* candidate consumes at most SampleLimit measurements and then strictly
\* narrows the safe/unsafe bracket.  A quiet view searches first for the
\* preferred steady cadence.  If quality debt remains and the independent
\* static deadline is available, it may advance exactly once to the static
\* goal, preserving the steady-safe lower bound.  A ready pose-only view may
\* start at the static goal: its retained presentation is the initial
\* candidate and may be coarsened only if it misses that hard deadline.  No
\* timer, repaint, or throughput-EMA update may reopen either search.
\*
\* This model proves ownership and termination.  It deliberately does not
\* model the numeric frame classifier or perceptual ordering of candidates;
\* those require independent property tests and graphical qualification.

EXTENDS FiniteSets, Naturals, TLC

CONSTANT CandidateCount, SampleLimit

ASSUME /\ CandidateCount > 0
       /\ SampleLimit > 0

Phases == {"choose", "measure", "goal_done", "terminal"}
Goals == {"steady", "static"}
Candidates == 1..CandidateCount

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
          certificateRevision,
          priorWidth,
          narrowed,
          goalTransitions

vars == <<phase, goal, trueSteadyCapacity, trueStaticCapacity, staticEligible,
          startAtStatic,
          safe, unsafe, candidate, samplesRemaining, measured,
          certificateRevision, priorWidth, narrowed, goalTransitions>>

Capacity(g) == IF g = "steady" THEN trueSteadyCapacity ELSE trueStaticCapacity

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
    /\ certificateRevision \in 0..1
    /\ priorWidth \in 1..(CandidateCount + 1)
    /\ narrowed \in BOOLEAN
    /\ goalTransitions \in 0..1

BracketSound == safe <= Capacity(goal) /\ Capacity(goal) < unsafe

CandidateOwned ==
    phase = "measure" =>
        /\ safe < candidate
        /\ candidate < unsafe
        /\ candidate \notin measured
        /\ samplesRemaining > 0

TerminalCertificate ==
    phase = "terminal" =>
        /\ unsafe = safe + 1
        /\ safe = Capacity(goal)
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
    /\ certificateRevision = 0
    /\ priorWidth = CandidateCount + 1
    /\ narrowed = FALSE
    /\ goalTransitions = IF startAtStatic THEN 1 ELSE 0

ChooseCandidate ==
    /\ phase = "choose"
    /\ unsafe > safe + 1
    /\ candidate' = safe + (unsafe - safe) \div 2
    /\ samplesRemaining' = SampleLimit
    /\ phase' = "measure"
    /\ narrowed' = FALSE
    /\ UNCHANGED <<goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, measured,
                    certificateRevision, priorWidth, goalTransitions>>

ConsumeSample ==
    /\ phase = "measure"
    /\ samplesRemaining > 1
    /\ samplesRemaining' = samplesRemaining - 1
    /\ narrowed' = FALSE
    /\ UNCHANGED <<phase, goal, trueSteadyCapacity, trueStaticCapacity,
                    staticEligible, startAtStatic, safe, unsafe, candidate, measured,
                    certificateRevision, priorWidth, goalTransitions>>

ClassifyCandidate ==
    /\ phase = "measure"
    /\ samplesRemaining = 1
    /\ priorWidth' = unsafe - safe
    /\ measured' = measured \cup {candidate}
    /\ IF candidate <= Capacity(goal)
          THEN /\ safe' = candidate
               /\ unsafe' = unsafe
          ELSE /\ safe' = safe
               /\ unsafe' = candidate
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
                    samplesRemaining, measured, priorWidth, goalTransitions>>

Next ==
    \/ ChooseCandidate
    \/ ConsumeSample
    \/ ClassifyCandidate
    \/ AdvanceToStaticGoal
    \/ PublishTerminalCertificate

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(ChooseCandidate)
    /\ WF_vars(ConsumeSample)
    /\ WF_vars(ClassifyCandidate)
    /\ WF_vars(AdvanceToStaticGoal)
    /\ WF_vars(PublishTerminalCertificate)

EventuallyCertified == <> (phase = "terminal")

=============================================================================
