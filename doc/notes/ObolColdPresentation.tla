---------------------- MODULE ObolColdPresentation --------------------------
\* Bounded cold-source presentation contract.
\*
\* This models the policy boundary for one large source.  A prompt local
\* source-order page is evidence that the worker is alive, but cannot replace
\* the whole-model overview: it is not spatially representative.  A global
\* preview may replace the overview only after it has complete coverage.  If
\* the source cannot obtain its bounded preparation reservation, the view
\* reaches an explicit overview-only constrained terminal state rather than a
\* blank frame or an unbounded retry loop.  Cache persistence is background
\* work after the view is ready and does not keep visual convergence pending.
\*
\* It deliberately abstracts mesh contents, wall-clock time, Qt scheduling,
\* and GPU work.  GUI evidence remains responsible for frame-time and image
\* quality; this model protects the discrete publication/liveness contract.

EXTENDS Naturals, TLC

SourceStates == {"new", "local", "global", "constrained", "complete"}
Presentations == {"overview", "mesh"}
Capacities == 0..1

VARIABLES capacity,
          source,
          presentation,
          visualPending,
          cachePending

vars == <<capacity, source, presentation, visualPending, cachePending>>

TypeOK ==
    /\ capacity \in Capacities
    /\ source \in SourceStates
    /\ presentation \in Presentations
    /\ visualPending \in BOOLEAN
    /\ cachePending \in BOOLEAN

Init ==
    /\ capacity \in Capacities
    /\ source = "new"
    /\ presentation = "overview"
    /\ visualPending = TRUE
    /\ cachePending = FALSE

PublishLocalEvidence ==
    /\ source = "new"
    /\ source' = "local"
    /\ UNCHANGED <<capacity, presentation, visualPending, cachePending>>

\* The local page is intentionally not a presentation replacement.
PublishGlobalPreview ==
    /\ source = "local"
    /\ capacity = 1
    /\ source' = "global"
    /\ presentation' = "mesh"
    /\ visualPending' = FALSE
    /\ cachePending' = TRUE
    /\ UNCHANGED capacity

\* Refusal is a terminal, diagnosable coverage state.  It cannot be retried
\* until an external capacity or view epoch creates a fresh source state.
Constrain ==
    /\ source = "local"
    /\ capacity = 0
    /\ source' = "constrained"
    /\ visualPending' = FALSE
    /\ UNCHANGED <<capacity, presentation, cachePending>>

CompleteCache ==
    /\ source = "global"
    /\ cachePending
    /\ source' = "complete"
    /\ cachePending' = FALSE
    /\ UNCHANGED <<capacity, presentation, visualPending>>

Next ==
    \/ PublishLocalEvidence
    \/ PublishGlobalPreview
    \/ Constrain
    \/ CompleteCache

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(PublishLocalEvidence)
    /\ WF_vars(PublishGlobalPreview)
    /\ WF_vars(Constrain)
    /\ WF_vars(CompleteCache)

\* The user always has a meaningful extent cue.  A local page cannot produce
\* the historical empty/fragmentary whole-object state.
NeverUndisplayable == presentation \in Presentations

LocalPageCannotHideCoverage ==
    source = "local" => presentation = "overview"

\* Mesh publication requires the globally representative transition, never
\* the local bootstrap alone.
MeshHasGlobalCoverage ==
    presentation = "mesh" => source \in {"global", "complete"}

\* A completed visual state is either a whole-object mesh or an explicit
\* capacity-constrained overview; background cache work is not visual debt.
TerminalVisualState ==
    ~visualPending =>
        (presentation = "mesh" \/ source = "constrained")

EventuallyVisualTerminal == <>(~visualPending)
EventuallyFullyQuiescent == <>(~visualPending /\ ~cachePending)

=============================================================================
