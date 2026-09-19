# Fairness and environment premises

Weak fairness in this suite is an implementation obligation, not a convenient
way to make TLC terminate.  A fair action is restricted to work BObol can
schedule after its guard remains continuously true.  User behavior, source
mutation, policy changes, view changes, cancellation, and input/session closure
are environment actions and are listed in `models.json`; the lint gate rejects
`WF_vars` applied to them.  It also rejects the former `WF_vars(Next)` pattern.

The production witnesses are:

| Fairness class | TLA+ action families | Concrete scheduling witness | Required failure edge |
|---|---|---|---|
| Owner-thread reducer | `Discover`, `Certify`, `Start*`, `Submit*`, `Apply*`, `Reject*`, `Plan*`, `Publish*`, `Retire*`, `DischargeOwner` | `BObolViewController::advanceProgressiveWork`, `BObolViewController::processPendingLodResults`, and the typed reducer owners named by each catalog record | A selected owner must expose an enabled typed transition; an owner/status bit alone is insufficient. |
| Worker completion | `Complete*`, `Finish*`, page/leaf construction, preparation units, compaction results | `BObolLodService` queues plus its `subscribeResultReady` callback, cache-writer completion, and owner-thread result drain | Potentially blocking external work needs a cancellation, timeout, or typed failure completion. Fairness must not assume an unbounded third-party call returns. |
| Frame acknowledgement | `CompleteFrame`, `CompleteCurrentRender`, `AcceptReport`, exact-presentation completion | The render endpoint's completed-frame report consumed by `BObolViewController`/`BObolWindowHost` | Deadline abort, stale report, endpoint loss, and successful current report are distinct outcomes and must all retire the in-flight owner. |
| Finite timer/deadline | debounce expiry, presentation deadline, compaction deadline, finite recovery | Host monotonic clock and the level-triggered pump request retained by `BObolWindowHost` and `BObolLodInteractionSession` | Timer expiry must retain a pump witness; clock arithmetic saturates and cannot silently wrap. |

Parameterized fair actions have the same witness for each finite member.  A
nondeterministic success/failure action is fair as a union: production promises
that the operation completes with a typed outcome, not that the environment
chooses success.  `ObolOccurrenceControl.CompleteLoad` is the reference form.

## Liveness closure rule

Liveness below an input boundary has one of these shapes:

- after an explicit final/closed revision, internal owners eventually reach a
  terminal outcome;
- after a quiet boundary, the epoch either settles or a later environment
  action interrupts it; or
- finite, monotonically bounded environment changes may occur, but no fairness
  forces them to occur.

An open interaction may therefore stutter forever without violating a
production liveness claim.  Once `CloseInput`, `EndInput`, `CloseCoverage`, or
the equivalent premise occurs, every remaining fair edge is mapped to one of
the concrete witnesses above.  Session teardown while a gesture is still open
is handled atomically by `ObolInteractionSession.CloseInput`; it is not hidden
behind fairness of a missing gesture-up event.

## Implementation audit rule

For every fair worker or frame action, the C++ audit must locate both the
normal completion and the path which fires when its executor, endpoint, or
consumer is torn down.  If neither a bounded timeout nor a cancellation/failure
completion exists, the corresponding TLA+ liveness result is conditional and
must not be reported as an implementation guarantee.
