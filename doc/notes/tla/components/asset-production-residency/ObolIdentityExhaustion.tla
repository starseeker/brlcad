----------------------- MODULE ObolIdentityExhaustion -----------------------
\* Machine-counter boundary for asynchronous authentication identities.
\*
\* The unbounded Naturals used by the surrounding models cannot establish
\* what a fixed-width C++ revision does at its maximum value.  This component
\* makes that boundary finite.  A requested semantic mutation either commits
\* with a fresh successor or enters a fail-stop terminal state; it never
\* reissues an earlier credential.  Diagnostic counters are outside this
\* contract because they do not authenticate work.

EXTENDS Naturals, TLC

CONSTANT MaxRevision

ASSUME MaxRevision \in Nat \ {0}

VARIABLES revision,
          evidenceRevision,
          mutationPending,
          accepted,
          halted,
          inputClosed,
          issuedRevisions

vars == <<revision, evidenceRevision, mutationPending, accepted, halted,
          inputClosed, issuedRevisions>>

TypeOK ==
    /\ revision \in 1..MaxRevision
    /\ evidenceRevision \in 0..MaxRevision
    /\ mutationPending \in BOOLEAN
    /\ accepted \in BOOLEAN
    /\ halted \in BOOLEAN
    /\ inputClosed \in BOOLEAN
    /\ issuedRevisions \subseteq 1..MaxRevision

Init ==
    /\ revision = 1
    /\ evidenceRevision = 0
    /\ mutationPending = FALSE
    /\ accepted = FALSE
    /\ halted = FALSE
    /\ inputClosed = FALSE
    /\ issuedRevisions = {1}

RequestMutation ==
    /\ ~inputClosed
    /\ ~halted
    /\ ~mutationPending
    /\ mutationPending' = TRUE
    /\ UNCHANGED <<revision, evidenceRevision, accepted, halted,
                    inputClosed, issuedRevisions>>

\* The semantic mutation and its fresh authentication identity commit in the
\* same transition.  MaxRevision is a valid final identity, not a wrap point.
Advance ==
    /\ mutationPending
    /\ revision < MaxRevision
    /\ revision' = revision + 1
    /\ mutationPending' = FALSE
    /\ accepted' = FALSE
    /\ issuedRevisions' = issuedRevisions \union {revision'}
    /\ UNCHANGED <<evidenceRevision, halted, inputClosed>>

\* Production terminates at this boundary.  The state-machine abstraction is
\* a closed, quiescent terminal state so TLC can check that the failed attempt
\* neither commits a mutation nor reuses a revision.
FailStop ==
    /\ mutationPending
    /\ revision = MaxRevision
    /\ mutationPending' = FALSE
    /\ accepted' = FALSE
    /\ halted' = TRUE
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<revision, evidenceRevision, issuedRevisions>>

CaptureEvidence ==
    /\ ~inputClosed
    /\ ~halted
    /\ evidenceRevision' = revision
    /\ accepted' = FALSE
    /\ UNCHANGED <<revision, mutationPending, halted, inputClosed,
                    issuedRevisions>>

AuthenticateEvidence ==
    /\ ~inputClosed
    /\ ~halted
    /\ evidenceRevision # 0
    /\ accepted' = (evidenceRevision = revision)
    /\ UNCHANGED <<revision, evidenceRevision, mutationPending, halted,
                    inputClosed, issuedRevisions>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ mutationPending' = FALSE
    /\ accepted' = FALSE
    /\ UNCHANGED <<revision, evidenceRevision, halted, issuedRevisions>>

Next ==
    \/ RequestMutation
    \/ Advance
    \/ FailStop
    \/ CaptureEvidence
    \/ AuthenticateEvidence
    \/ CloseInput

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(Advance)
    /\ WF_vars(FailStop)
    /\ WF_vars(AuthenticateEvidence)

AcceptedEvidenceIsCurrent ==
    accepted => evidenceRevision = revision

NoIdentityReuse ==
    issuedRevisions = 1..revision

MutationHasProgressWitness ==
    mutationPending => ENABLED Advance \/ ENABLED FailStop

HaltedIsQuiescent ==
    halted => inputClosed /\ ~mutationPending /\ ~accepted

HaltedOnlyAtExhaustion ==
    halted => revision = MaxRevision

Terminal == inputClosed /\ ~mutationPending

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyResolveMutation ==
    [](mutationPending => <>~mutationPending)

=============================================================================
