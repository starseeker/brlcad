# BObol formal terminology

This vocabulary applies to TLA+ models, architecture documents, control
diagnostics, and executable refinement tests.  A focused model may use a
smaller abstraction, but must document any intentionally narrower meaning.

## State and identity

- **Demand:** the desired semantic result supplied by a user, database,
  policy, or downstream owner.  Demand is a value, not merely a wakeup.
- **Revision:** a monotonically advancing identifier for one semantic domain.
  It changes exactly when that domain's meaning changes.
- **Epoch:** the interval during which a selected set of revisions is frozen
  for one transaction.  It is not a synonym for every individual revision.
- **Evidence stamp:** the complete typed revision tuple consumed by a plan or
  certificate.  BObol's canonical planning stamp contains inventory,
  availability, visibility, view, policy, and capacity revisions.
- **Source-route revision:** the lifetime of the provider/cache path
  authorized to produce an asynchronous result.  It is independent of both
  the semantic population epoch and the current demand revision.
- **Current:** certified for the complete evidence stamp now owned by the
  consumer.
- **Stale:** no longer certified for the current evidence stamp.
- **Superseded:** retired because a newer semantic demand owns the boundary.
  Supersession is not a provider failure.
- **Cancelled:** explicitly prevented from completing or publishing.  A
  cancelled worker result may still need mechanical retirement.

## Work and ownership

- **Obligation:** a standing semantic fact which must be discharged before a
  foreground terminal outcome.
- **Owner:** the one protocol component authorized to choose the obligation's
  successor transition.
- **Progress witness:** a worker, cursor, owner-thread pump, frame
  acknowledgement, or finite timer which can actually advance an obligation.
  An owner label or notification edge alone is not a witness.
- **Work unit:** one member of a finite mechanical rank.  Consuming it changes
  progress without independently choosing policy.
- **Candidate:** an uncommitted possible successor.
- **Plan:** a deterministic, revision-stamped policy decision and its complete
  successor evidence.
- **Cursor:** bounded resumable mechanical progress through an accepted plan;
  it contains no independent quality or terminal policy.
- **Result:** immutable producer output submitted to its consumer for current-
  demand validation.
- **Certificate:** revision-bound evidence authorizing reuse, constraint, or a
  terminal conclusion.
- **Transaction:** the sole atomic owner of a multi-step commit or handoff.

## Representation and publication

- **Coverage:** a truthful representation exists for every required visible
  occurrence, possibly using a structural fallback.
- **Availability:** immutable producer data exists and can be considered by
  planning; it is not necessarily resident or displayed.
- **Residency:** the admitted data needed for a selected representation is in
  the live working set.
- **Retained-scene commit:** validated candidate data atomically becomes the
  renderer-visible retained scene.
- **Presentation commit:** an accepted allocation or mutation becomes the
  target of an exact framebuffer transaction.
- **Frame acceptance:** a completed report for the exact current target
  becomes the accepted framebuffer evidence.
- **Durable-cache commit:** complete immutable producer data becomes reusable
  across process lifetimes.

Use the qualified terms above instead of bare **publication** whenever more
than one boundary could be meant.

## Outcomes and verification

- **Ready:** the requested foreground contract is satisfied and no foreground
  obligation remains.
- **Constrained:** unmet quality demand is paired with explicit current
  evidence explaining why no admissible successor exists.
- **Failed:** the current demand cannot produce a valid result.  Failure is not
  interchangeable with cancellation or supersession.
- **Terminal:** ready, constrained, or failed with no foreground owner.  Some
  subordinate background reclamation may remain when explicitly modeled.
- **Safety invariant:** a state predicate which must hold in every reachable
  state.
- **Liveness property:** a temporal property asserting eventual progress under
  documented input-closure and fairness assumptions.
- **Fairness assumption:** an explicit scheduler/environment premise.  It may
  not assume away a production wakeup or progress obligation.
- **Refinement mapping:** the documented and executable projection from
  concrete C++ facts/transitions to an abstract model.
- **Bounded model check:** TLC's exhaustive exploration of one configured
  finite state space.  It is evidence for the modeled abstraction, not a proof
  of unbounded implementation correctness.
