----------------------- MODULE ObolCadTimingEvidence -----------------------
\* A retained CAD work record and its duration are one threshold-stamped
\* sample.  Non-CAD host frames may update unrelated timing, but cannot drive
\* the finite point-quality bracket.  MaxCut + 1 represents a renderer whose
\* irreducible occurrence population misses the preferred cadence.

EXTENDS Integers, Naturals

CONSTANT MaxCut

Cuts == 1..MaxCut
Outcomes == {"none", "safe", "overload"}

VARIABLES cut, unsafeCut, safeCut, sampleCut, sampleOutcome,
          hostOutcome, sustainableCut, initialFramePending,
          replayPending, terminal

vars == <<cut, unsafeCut, safeCut, sampleCut, sampleOutcome,
          hostOutcome, sustainableCut, initialFramePending,
          replayPending, terminal>>

TypeOK ==
    /\ MaxCut > 1
    /\ cut \in Cuts
    /\ unsafeCut \in 0..MaxCut
    /\ safeCut \in 0..MaxCut
    /\ sampleCut \in 0..MaxCut
    /\ sampleOutcome \in Outcomes
    /\ hostOutcome \in {"safe", "overload"}
    /\ sustainableCut \in 1..(MaxCut + 1)
    /\ initialFramePending \in BOOLEAN
    /\ replayPending \in BOOLEAN
    /\ terminal \in BOOLEAN

Init ==
    /\ cut = 1
    /\ unsafeCut = 0
    /\ safeCut = 0
    /\ sampleCut = 0
    /\ sampleOutcome = "none"
    /\ hostOutcome = "safe"
    /\ sustainableCut \in 1..(MaxCut + 1)
    /\ initialFramePending = TRUE
    /\ replayPending = FALSE
    /\ terminal = FALSE

\* The first exact retained framebuffer proves presentation correctness, but
\* may still contain one-time upload or command preparation.  An event-driven
\* host must turn that missing timing evidence into one explicit replay owner;
\* an unrelated repaint is not an implicit successor.
CompleteInitialCadFrame ==
    /\ ~terminal
    /\ initialFramePending
    /\ sampleCut = 0
    /\ initialFramePending' = FALSE
    /\ replayPending' = TRUE
    /\ UNCHANGED <<cut, unsafeCut, safeCut, sampleCut, sampleOutcome,
                    hostOutcome, sustainableCut, terminal>>

\* The unchanged replay publishes work identity, duration class, and
\* classifier as one record.  The controller may consume it exactly once.
CompleteReusableCadReplay ==
    /\ ~terminal
    /\ ~initialFramePending
    /\ replayPending
    /\ sampleCut = 0
    /\ sampleCut' = cut
    /\ sampleOutcome' = IF cut >= sustainableCut
                         THEN "safe" ELSE "overload"
    /\ replayPending' = FALSE
    /\ UNCHANGED <<cut, unsafeCut, safeCut, hostOutcome,
                    sustainableCut, initialFramePending, terminal>>

\* Faceplate, overlay, and cached-image frames deliberately do not mutate the
\* CAD sample, even when their duration class differs from it.
RecordNonCadFrame ==
    /\ ~terminal
    /\ hostOutcome' = IF hostOutcome = "safe"
                      THEN "overload" ELSE "safe"
    /\ UNCHANGED <<cut, unsafeCut, safeCut, sampleCut, sampleOutcome,
                    sustainableCut, initialFramePending, replayPending,
                    terminal>>

ConsumeOverload ==
    /\ ~terminal
    /\ sampleCut = cut
    /\ sampleOutcome = "overload"
    /\ LET nextUnsafe == IF cut > unsafeCut THEN cut ELSE unsafeCut
           adjacentSafe == safeCut > nextUnsafe /\
                           safeCut - nextUnsafe <= 1
       IN
       /\ unsafeCut' = nextUnsafe
       /\ sampleCut' = 0
       /\ sampleOutcome' = "none"
       /\ IF cut = MaxCut
             THEN /\ cut' = cut
                  /\ replayPending' = FALSE
                  /\ terminal' = TRUE
          ELSE IF adjacentSafe
             THEN /\ cut' = safeCut
                  /\ replayPending' = FALSE
                  /\ terminal' = TRUE
          ELSE /\ cut' \in IF safeCut > nextUnsafe
                              THEN (nextUnsafe + 1)..safeCut
                              ELSE (cut + 1)..MaxCut
               /\ replayPending' = TRUE
               /\ terminal' = FALSE
    /\ UNCHANGED <<safeCut, hostOutcome, sustainableCut,
                    initialFramePending>>

ConsumeSafe ==
    /\ ~terminal
    /\ sampleCut = cut
    /\ sampleOutcome = "safe"
    /\ LET nextSafe == IF safeCut = 0 \/ cut < safeCut THEN cut ELSE safeCut
           adjacentUnsafe == unsafeCut > 0 /\ nextSafe > unsafeCut /\
                             nextSafe - unsafeCut <= 1
       IN
       /\ safeCut' = nextSafe
       /\ sampleCut' = 0
       /\ sampleOutcome' = "none"
       /\ IF cut = 1 \/ adjacentUnsafe
             THEN /\ cut' = cut
                  /\ replayPending' = FALSE
                  /\ terminal' = TRUE
          ELSE /\ cut' \in IF unsafeCut > 0
                              THEN {unsafeCut +
                                  (nextSafe - unsafeCut) \div 2}
                              ELSE 1..(cut - 1)
               /\ replayPending' = TRUE
               /\ terminal' = FALSE
    /\ UNCHANGED <<unsafeCut, hostOutcome, sustainableCut,
                    initialFramePending>>

Next ==
    \/ CompleteInitialCadFrame
    \/ CompleteReusableCadReplay
    \/ RecordNonCadFrame
    \/ ConsumeOverload
    \/ ConsumeSafe

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(CompleteInitialCadFrame)
    /\ WF_vars(CompleteReusableCadReplay)
    /\ WF_vars(ConsumeOverload)
    /\ WF_vars(ConsumeSafe)

SampleIsCoherent ==
    /\ (sampleCut = 0) = (sampleOutcome = "none")
    /\ (sampleCut > 0 =>
           sampleOutcome = IF sampleCut >= sustainableCut
                           THEN "safe" ELSE "overload")

BracketIsOrdered ==
    /\ (unsafeCut > 0 => unsafeCut < sustainableCut)
    /\ (safeCut > 0 => safeCut >= sustainableCut)
    /\ (unsafeCut > 0 /\ safeCut > 0 => unsafeCut < safeCut)

TimingOwnerIsExclusive ==
    /\ ~(initialFramePending /\ replayPending)
    /\ ~(initialFramePending /\ sampleCut > 0)
    /\ ~(replayPending /\ sampleCut > 0)

TerminalHasNoSample == terminal =>
    /\ sampleCut = 0
    /\ ~initialFramePending
    /\ ~replayPending

ActiveHasTimingWitness == ~terminal =>
    initialFramePending \/ replayPending \/ sampleCut > 0

EventuallyTerminal == <>terminal

=============================================================================
