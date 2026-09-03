# libBObol TLA+ to C++ conformance audit

Last reviewed: 2026-09-03

This is the implementation-facing conformance ledger for the BObol formal
suite.  `tla/models.json` remains the authoritative proof graph and test
catalog; this document records whether the production owners actually refine
the modeled transitions.  A TLC pass proves the bounded abstract relation,
not the C++ mapping.  A row is conforming only when the formula, production
owner, and executable witness agree.

The audit began with the failure families in
`libbobol_engineering_lessons.md`, then considered adjacent races rather than
limiting its scope to failures already observed.  In particular, address and
dense-slot reuse motivated independent identity domains; stale demand failures
motivated authentication before outcome interpretation; lost wakeups motivated
level-triggered owners; and mutable shared generations motivated explicit
consumer lifetime rather than reliance on pointer ownership.

## Status vocabulary

- **Conforming**: the model-to-C++ mapping is explicit and an executable test
  exercises the refinement boundary.
- **Partial**: the core owner exists, but an identified transition, trace, or
  failure edge lacks executable evidence.
- **Nonconforming**: the current C++ ownership permits a behavior forbidden by
  the model.
- **Outside abstraction**: retain a separate data-plane, fault-injection,
  sanitizer, performance, image, or platform gate.

## Shared terminology

| Formal term | Production meaning | Authoritative C++ representation |
|---|---|---|
| evidence stamp | all semantic inputs which authorize one admission plan | `BObolLodAdmissionRevisionStamp`: inventory, availability, visibility, view, policy, capacity |
| demand | the view-local request which authorizes presentation | `BObolViewEpoch` plus `BObolPolicyEpoch`; neither is asset identity |
| source route | lifetime of the source owner to which a result returns | `BObolSourceRoutingId` |
| source population | lifetime of the compact dense registry addressed by an entry index | `BObolSourcePopulationEpoch` |
| current result | route, population, view, and policy all match their current owners | `BObolLodResultAuthenticationContract` |
| publish | install a useful current result | `BObolLodResultDisposition::PUBLISH` followed by `BObolViewLodState::consumeSourceResult` |
| terminal failure | exact-current provider/cache failure evidence | `RECORD_TERMINAL_FAILURE` and the demand-scoped CAD occurrence failure map |
| retry | no failure evidence; recreate current demand | `RETRY_CURRENT_DEMAND` and the controller demand replay edge |
| supersede | discard an identity-mismatched result before interpreting its outcome | `SUPERSEDE`; no payload or failure-map mutation |
| work owner | the one enabled finite mechanism which can discharge an obligation | `BObolLodControlRefinement::Owner` |
| consumer lease | one live view/generation's claim on an in-flight shared producer or queued result | `BObolSharedProducer::leases`, keyed by service generation and carrying that consumer's newest demand |
| authentication identity | a fixed-width credential whose equality authorizes work, publication, or ownership | checked strong values and `identity_counter_private.h`; exhaustion is fail-stop and never wraps |
| diagnostic counter | an observation total which authorizes no behavior | saturating counter; it may stop increasing but may not be used as an identity |

An occurrence-targeted result must carry nonzero route and population
credentials even when a focused test uses semantic-key routing with
`sourceEntryIndex == UINT32_MAX`.  Empty credentials are not a wildcard.  A
source-wide result has an empty occurrence key and is outside compact
occurrence authentication.

## Conformance by proof boundary

| Proof boundary | Production enforcement | Executable evidence | Status |
|---|---|---|---|
| `ObolProgressivePipeline`: one six-domain stamp, one bounded cursor, one finite owner, derived terminal outcome | `lod_revision_private.h`, `lod_control_private.h`, coordinator evidence values and controller effect executors; exact Obol GPU keys; exact/interned capacity, structural-frontier, and compact-presentation credentials | `test_lod_coordinator.cpp` six-domain revision, owner, witness, typed-outcome, equivalent-population, and collision tests; retained-allocation identity tests; randomized completed-pass and transition-journal tests; focused Obol renderer tests | **Conforming** after `CXX-STAMP-001`: partial and caller-assembled planning stamps are unrepresentable, every domain independently invalidates stored plans and allocation keys, and mutable-state cache authorization retains exact inputs rather than a folded digest. |
| `ObolControlRefinement`: 29 independent facts refine to nine work classes and one fixed-precedence owner | `BObolLodControlRefinement`, `BObolLodControlTransitionScope`, controller transition journal, and runtime violation projection | per-field and exhaustive owner combinations in `test_lod_coordinator.cpp`; randomized journal test; complete graphical trace checker | **Conforming** after `CXX-TRACE-001`: every retained transition has typed endpoints and the audit sentinel rejects an unmapped writer. |
| `ObolResultAuthentication`: independent route, population, and demand authentication before publish/failure/retry | service result normalization, `lod_result_authentication_private.h`, controller publication gate, view-state reducer | direct disposition table plus the 50-case asynchronous controller matrix in `test_lod_update_action.cpp`; service request-matching tests | **Conforming** after `CXX-RP-001`. |
| `ObolIdentityExhaustion`: a fixed-width identity advances without reuse or fail-stops before a semantic mutation | `BObolLodStrongUInt64`, libBObol `identity_counter_private.h`, and Obol `CadIdentityCounter.h` across retained-plan, classifier, preparation, timing, and resource evidence | successor boundary tests in `test_lod_coordinator.cpp` and Obol `test_cad_instance_records.cpp`; cross-session trace-token regression in `test_window_host.cpp` | **Conforming** after `CXX-ABA-001`. |
| `ObolSharedAssetLease`: a surviving or late consumer retains shared build/result ownership; only the final lease cancels | typed shared-producer state, generation leases, producer-aware cancellation, and per-generation payload/replay drains in `BObolLodService` | same-generation coalescing plus six cross-generation build/result/cancellation lifecycle cases in `test_lod_service.cpp` | **Conforming** after `CXX-SL-001`. |
| admission/allocation/capacity and terminal-quality component models | typed coordinator policies and retained-allocation transaction | focused exhaustive/property tests and renderer evidence tests | **Conforming** for modeled control relations; numeric/visual sufficiency remains outside TLA+. |
| retained CAD mutation and exact-frame compositions | private staging, checked Obol mutation, exact target/report revisions | Obol mutation tests, CAD presentation/frame tests, controller endpoint and commit-denial tests | **Conforming** for atomic precommit resource denial after `CXX-TXN-001` and typed endpoint cancellation after `CXX-LIFE-001`. |
| host, policy, interaction, and semantic-presentation lifecycle | typed host request, interaction session, exact-frame and policy retirement ledgers | off/on/off, interrupted frame, callback/unsubscribe, endpoint-loss, and semantic-only presentation tests | **Conforming** after `CXX-LIFE-001`: endpoint loss atomically retires controller and host owners, while callback removal drains or defers every live dispatch. |
| cache publication, durability, geometry truth, renderer pixels, resource magnitude | staged payload records, atomic name mapping, provider validation, immutable geometry, renderer implementations | denied replacement, corrupt-cache, geometry, image, GUI, sanitizer, and performance matrices | **Conforming** for modeled atomic name publication; crash consistency, numeric truth, and resource magnitude remain outside the abstraction as noted in `tla/RISK_COVERAGE.md`. |

## Closed vertical finding: result publication

### CXX-RP-001 — stale result could bypass complete authentication

Before this pass, compact route and population were checked separately in the
view state, while view/policy arbitration in the controller allowed an old
mesh to publish as a cold bootstrap when no equal or newer mesh was resident.
That violated `ObolResultAuthentication.StaleResultCannotApply`.  It also
interpreted some provider outcomes before proving identity, so an obsolete
terminal error could become failure evidence for a replacement occurrence.

The production boundary now performs one allocation-free decision over three
independent identity domains.  Any mismatch yields `SUPERSEDE` before provider
status is interpreted, causes no payload or failure-map mutation, and arms one
current-demand replay.  Current `READY` publishes; current `CACHE_MISS`,
`STALE`, and `ERROR` record exact-demand terminal failure; other provider
statuses retry without failure evidence.

The cumulative progressive-mesh rebase remains legal only as a
pre-authentication refinement.  It must prove identical immutable asset
ownership and reconstruct the current page set, counts, and prepared renderer
generation.  Only then may it rewrite the request to current demand and enter
the publisher.  It is not permission to apply a stale result.

`ObolResultAuthentication` now includes the same publish, terminal-failure,
and retry dispositions.  TLC checks that a stale result enables neither
publication nor failure recording and cannot create an accepted outcome.  The
configured finite run explores one bounded retry because production retry
count is deliberately not a formal resource claim.

The C++ table tests every provider status.  The asynchronous matrix crosses
all ten statuses with current identity and independently stale route,
population, view, and policy identities.  Stale rows run while the occurrence
is cold, so a future bootstrap exception cannot silently restore the defect.
The stricter contract also corrected hand-built fixtures which had omitted
route/population credentials that production requests always carry.

## Closed vertical finding: shared producer lifetime

### CXX-SL-001 — coalescing is not a consumer lease

Before this pass, `lod_request_active_key` correctly removed occurrence identity when
`coalesceAssetProducer` is true, and it preserves route and compact-population
identity.  `activeRequestKeyCounts` and `latestActiveRequests` then serialize
the producer and retain its newest demand.  This is request coalescing, not the
multi-consumer lifetime required by `ObolSharedAssetLease`:

- a duplicate submission returns no task ID and records no consumer or
  generation lease;
- the surviving work item and completed result retain only the first task's
  generation;
- `cancelGeneration` cancels that work and erases that generation's queued
  result without consulting another generation which requested the same
  asset; and
- `drainGenerationResults` can deliver the result only to its original
  generation.

Consequently, closing the first view could cancel a build still needed by a
second view, and a late join has no durable claim on a queued result.  The
post-publication `residentMeshConsumerDemands` snapshots protect compaction;
they do not lease an active producer and cannot satisfy this model.

The service now owns one typed producer record per stable active asset key and
one lease per service generation.  Duplicate discovery in both submission and
the earlier demand-update fast path records the lease and that consumer's
newest demand.  Cancelling a generation removes only its lease; worker,
preview, debug-delay, and cache cancellation consult the producer lifetime and
continue while any lease survives.  Pending work and undelivered results
retire only with the final lease.

Only the producer generation receives its demand-specific payload.  Other
generations receive a lightweight, exact-current `SUPERSEDED` replay result,
which the result-authentication reducer maps to `RETRY_CURRENT_DEMAND`.  This
avoids copying large mesh arrays or publishing geometry prepared under another
view's admission certificate; the awakened consumer instead binds the shared
resident asset under its own current demand.  Late joins work during both the
build and queued-result states.

Six executable lifecycle cases cover two generations sharing one delayed
build, first-generation cancellation while the second survives, late join
during build, late join while a result is queued, final-consumer cancellation
during build, and final-consumer cancellation with a queued result.  They
assert one provider call, generation-scoped delivery or replay, and complete
retirement of producer, lease, request, result, and in-flight ownership.
`ObolSharedAssetLease.ResultRetainedForUnsatisfiedConsumer` makes the queued
result rule explicit.  The affected TLC closure passed with baseline counts:
237 generated / 96 distinct lease states, 7,067 / 1,900 asset-publication
states, 712,617 / 356,352 control-refinement states, and 10,167,666 /
4,303,828 canonical pipeline states.

## Closed vertical finding: complete control-transition refinement

### CXX-TRACE-001 — sampled states were not an execution trace

The prior qged observer sampled distinct convergence states at event-loop
intervals.  Its projection and cycle checks were useful, but several reducer
effects could occur between observations.  A passing report therefore proved
only that the states it happened to see were valid; it could not prove that
each concrete transition belonged to the finite refinement relation.

The controller now owns an opt-in bounded journal.  Every retained record has
an immutable before and after endpoint, a contiguous serial, and one event
from the external-input plus nine-owner alphabet.  Owner-thread effect
executors use nested `BObolLodControlTransitionScope` boundaries.  Cohesive
external-input setters suppress their implementation-detail nested scopes so
revision, obligation, and successor wake publication remain one atomic
transition.  Worker result callbacks publish only a typed publication event;
the owner captures its endpoints without reading Coin state from the worker.

The endpoint includes the complete control fact/obligation/owner projection,
all admission revisions, plan and presentation identities, host work, and the
exact renderer-preparation target signature and remaining-unit rank.  When
observable state changes between registered scopes, the journal emits
`unnamed`; it never guesses an owner.  Full buffers report drops rather than
overwriting history.  Neither the journal nor its checker is consulted by
production policy.

The qged checker now rejects an unnamed or unknown event, noncontiguous
serial, missing endpoint field, discontinuous before/after pair, event/owner
mismatch, rank regression, local invariant, and trace truncation.  Its
adversarial shell fixtures cover each new rejection.  A 512-step deterministic
controller trace checks boundary nesting, continuity, and bounded-loss
reporting in C++.

The first complete OSMesa draw/settle/zoom trace exposed a real missing
refinement witness: consuming a queued render cleared `RENDER` before the
exact-frame reducer completed.  `ObolHostWork` already models that interval as
`renderInFlight`; production now publishes the matching typed claimed-frame
host level and presentation witness and retires it only at completed or
interrupted-frame reduction.  The repaired representative trace contains 57
transitions, 15 claimed-frame endpoints, no unnamed events or drops, and
passes the independent refinement checker.

## Closed vertical finding: complete evidence stamps

### CXX-STAMP-001 — admission stamps permitted partial construction

The canonical TLA+ model already represents both current and in-progress
planning identity as total six-field records.  C++ did not enforce that
shape: `BObolLodAdmissionRevisionStamp` was a default-constructible aggregate,
and its coordinator snapshot and tests assigned companion fields one at a
time.  A future asynchronous constructor could therefore omit visibility or
capture fields on opposite sides of a semantic mutation while still passing
the type checker.

The stamp is now an immutable, allocation-free value requiring typed
inventory, availability, visibility, view, policy, and capacity arguments.
The only empty value is the named `administrative()` sentinel used by
unstamped immediate reducer results and inactive certificates.  Stored plan,
cursor, allocation, constraint, and acceptance values initialize that
sentinel explicitly; asynchronous inputs and results can replace it only by
copying a complete stamp.  The owner-thread coordinator is the sole production
constructor from live revision owners and captures all six in one expression.

All six mutation domains now have executable stale-identity checks for the
admission plan, its resumable cursor, and the retained-allocation input key.
Compile-time checks reject a default stamp and a five-domain constructor while
preserving trivial-copy and fixed-size behavior.  The audit also found three
intentional partial comparisons for renderer-capacity evidence.  They now use
one named `sameRendererCapacityProblem` projection: inventory, availability,
view, and policy define that measured problem; visibility is a replaceable
allocation input and capacity is the measurement's output.  This distinction
prevents an ad hoc partial comparison from masquerading as full plan
authentication.

No formula change was needed: `ObolProgressivePipeline` already requires the
complete stamp at plan start, each planning step, completion, and stale abort.
The affected coordinator, retained-allocation, host, and trace tests pass.

The cross-repository audit subsequently found lower-level Obol GPU caches which
folded several revisions, progressive ranges, or flattened-atlas credentials
into a 64-bit digest and then treated equality as proof that an old upload was
current.  That conflicts with the model even though collision is unlikely:
probability is not an authorization invariant.  Aggregate-proxy uploads now
retain the exact source, attribute, presentation, draw-mode, and software-path
tuple.  Progressive cut buffers retain their exact representation, lineage,
interval, quantization, and packed range vector.  Flattened wire and shaded
atlas entries retain instance, part, cut, source revision, transform, and
source interval.  Hashing remains legal only as a lookup accelerator followed
by exact equality.

The same rule applies above the renderer.  Retained-allocation capacity reuse
now carries a non-reused population identity issued only after the view owner
compares the exact ordered selected population.  Equal diagnostic digests do
not authorize reuse.  Structural-frontier deduplication uses an exact bounded
interner; eviction can cause only a redundant rebuild, never false reuse.
Compact-source and cross-source presentation caches compare exact revision and
membership vectors, including semantic revision, rather than structure/style
digests.  Length-prefixed layered payload keys and exact stable-handle
interning close the corresponding variable-string identity paths.

Completed-frame aggregation follows the same rule.  Obol assigns every
`SoCADAssembly` a process-lifetime identity, including across address reuse.
libBObol canonicalizes the exact `(assembly identity, local serial)` set and
issues a checked non-reused token for draw execution, GPU timing, and GPU
resource observations.  Invalid or absent evidence retires the current set, so
an `A -> absent -> A` source sequence cannot reuse its earlier aggregate token.
Obol likewise replaces the saturating sum of its assembly- and renderer-owned
preparation serials with an exact-pair change identity.

No TLA+ formula changed for these closures.  The formulas already compare
abstract records and identities exactly; the C++ digest shortcuts were invalid
refinements of that relation.  Strong 128-bit immutable geometry identifiers
and persistent mesh content hashes remain deliberate content-address domains,
not mutable evidence stamps.  Their byte-level validity and cache durability
remain executable/operational obligations rather than claims of the bounded
control models.

## Closed vertical finding: transactional publication

### CXX-TXN-001 — denial could notify or erase before commit

The retained-scene path had pure staging and checked mutation calls, but an
outer `SoCADAssembly` update scope was opened before the final resource gate.
An injected denial changed no scene bytes yet still closed that empty scope and
advanced the node notification identity.  Resource admission now occurs after
all compact-presentation preflight and before the first update scope.  The
executable fixture proves that both retained-scene and presentation denial
preserve the prior scene and `SoNode` notification ID, while the succeeding
commit changes both together.

The durable mesh-LoD path had the inverse problem.  Forced refresh deleted the
old name-to-content mapping before candidate generation.  A write, allocation,
or disk failure could therefore leave a complete immutable predecessor on disk
but make it undiscoverable.  Refresh and authored-mesh replacement now build
all content-addressed records first and swap the name mapping only after the
candidate is complete.  Deterministic failure points abort component and
spatial transactions immediately before commit and reject the final mapping
write.  Tests cover both no-prior-object denial and failed replacement of a
known-good hierarchy; the latter must retain the same key and reopenable
payload.  Optional in-process name hints are cleared on post-commit allocation
failure so they cannot mask the authoritative disk mapping.

This finding required a formula correction.  The former
`ObolLiveSpatialPublication` and asset composition represented durability with
only a Boolean marker, so they could not express replacement of an existing
mapping—the exact state lost by the C++ bug.  Both models now carry a baseline
mapping, an atomic candidate publication, and `DenyCacheCommit`.
`DeniedCacheCommitPreservesMapping` proves that denial leaves the baseline
discoverable, while `CacheMarkerRequiresFinalHierarchy` still prevents early
publication.  TLC checks the denial action as required coverage in both the
focused model and its parent composition.

The guarantee is deliberately at the precommit boundary.  Complete-scene
replacement uses a staged no-throw swap.  Sparse retained mutation reserves
for its bounded journal but does not clone a 150k scene to recover from
process-wide allocator exhaustion during the mechanical commit.  Filesystem
crash consistency and kill-at-instruction durability remain separate
operational tests rather than claims of this bounded control model.

## Closed vertical finding: endpoint and callback lifetime

### CXX-LIFE-001 — endpoint loss left work without a host

The display endpoint previously detached its callback and host without
retiring work retained by a borrowed controller.  An open gesture could remain
in its interaction/debounce protocol, a consumed render request could remain a
claimed-frame witness, and worker, coverage, presentation, or capacity owners
could remain live even though no endpoint existed to advance them.  This was
not an orderly producer close: the remaining frame and timer fairness premises
were false after the host disappeared.

Endpoint loss is now one explicit cancellation transaction.  Before clearing
the host levels it resets renderer-derived presentation controls, cancels the
active service generation, retires every automatic-LoD domain, closes the
interaction session, invalidates capacity evidence, and clears pending and
claimed frames plus the shared pump.  The postcondition asserts both the
complete automatic-domain retirement predicate and an empty typed host-work
snapshot.  Immutable resident payloads and user policy survive; binding a new
endpoint starts a fresh coverage/capacity transaction and reprojects any
independent provider work.

The same audit found that presentation callback removal copied a raw callback
and user pointer under a mutex, then invoked it without a dispatch lease.
Removal could return while the callback still used freed endpoint state.
Frame callbacks had a dispatch count, but self-unsubscribe waited for its own
dispatch and deadlocked.  Both callback kinds now use the same reference-
counted close protocol.  External removal waits for all dispatches; removal
from the current or any enclosing nested callback defers deletion until that
dispatch unwinds.  An allocation failure preserves the prior live callback
rather than silently dropping it.

Executable cases close an endpoint around a borrowed controller with an open
gesture, active service generation, progressive work, and claimed capacity
frame, and require zero surviving host flags, control facts, obligations,
owner, generation, or capacity search.  Separate cases prove external
presentation teardown blocks until callback return, self-unsubscribe is safe,
and a nested callback can unsubscribe its outer dispatch without deadlock.
The pre-existing factory race still proves host destruction waits for an
active frame callback.

`ObolInteractionSession.CloseInput` already modeled close during a gesture.
`ObolHostWork` now distinguishes endpoint loss from drainable
`CloseProducers`: `CloseEndpoint` atomically cancels notification, loop, pump,
pending render, claimed render, exact-frame, publication-timer, source, and
provider ownership.  It is required nonzero action coverage and exports
`EndpointClosureRetiresAllWork` to the lifecycle composition contract.  The
focused TLC run passes at 89,521 generated / 12,924 distinct states to depth
17.

## Closed vertical finding: fixed-width identity exhaustion

### CXX-ABA-001 — counter wrap could authenticate retired evidence

The formal revision domains were monotonic naturals, but several production
counters were fixed-width integers.  Most skipped zero after unsigned wrap;
`compact_next_revision` explicitly returned to one.  A sufficiently old
revision, generation, route, handle, lineage, plan serial, or callback token
could therefore compare equal to a new credential and authenticate retired
work.  Skipping the zero sentinel reduced the collision set by one value; it
did not prevent ABA reuse.

All increment sites which create authentication, ordering, or ownership
credentials now use the same checked-successor contract in their owning
repository.  It covers libBObol local and atomic counters, zero-valid
sequences, nonzero identities, public strong LoD epoch values, and Obol's CAD
plan, classifier, preparation, GPU-resource, and timing evidence.  A current
identity may advance to the maximum value; a stored next-identity allocator
reserves that value as its exhaustion marker.  In either representation, an
operation which needs an unavailable successor terminates before it can
commit the associated semantic mutation.  This fail-stop boundary is
intentional: most affected APIs cannot roll back every already-observed
dependent object, and continuing after exhaustion would be less safe than
stopping.  Obol resource teardown no longer resets externally observed upload
or timer allocation domains, so context recreation cannot manufacture an ABA
without numeric overflow.  Diagnostic totals are explicitly separate and
saturate because they authorize no behavior; resource byte/count arithmetic
remains governed by its own saturation and reservation tests.

Aggregate identities are checked at the same boundary.  An assembly object has
a process-lifetime identity rather than relying on its reusable address, and
multi-assembly execution, GPU timing, and resource samples are interned from
their exact canonical source tuples.  Obol's combined preparation token is
also derived from exact component equality rather than a saturating sum.

The audit also found an ABA path without numeric overflow.  Disabling and
re-enabling controller transition tracing reset its token allocator while an
old RAII scope could still exist.  A late destructor could then carry the same
token as a new scope and close the wrong transition frame.  Trace records may
be cleared between sessions, but their token and serial domains now remain
monotonic for the controller lifetime.  A regression leaves an old scope
alive across re-enable and proves that it is rejected rather than accepted as
the current scope.

`ObolIdentityExhaustion` closes the TLA+ abstraction gap with a finite revision
domain.  A requested mutation either atomically publishes a fresh successor
or takes `FailStop` at the boundary; issued revisions are never reused, stale
evidence cannot authenticate, and the halted state is closed and quiescent.
The focused TLC run explores the complete configured graph at 91 generated /
37 distinct states to depth 8.  Pure C++ boundary tests in both repositories
cover zero, the final valid successor, exhausted detection, local/atomic
allocation, and diagnostic saturation without deliberately terminating the
test process.

## Remaining audit findings

All identified formal-to-C++ vertical findings are closed.  The independent
numeric, geometry, image, concurrency, durability, performance, and platform
qualification gates below remain necessary; they are not downgraded by this
ledger status.

## Change gate

Any change to result routing, shared producer lifetime, revision identity, or
terminal outcome must update together:

1. the focused component model and applicable composition;
2. its production typed value/reducer rather than a parallel latch;
3. the executable table or lifecycle matrix at the real asynchronous seam;
4. `tla/models.json` ownership/test metadata and this ledger; and
5. the focused TLC baseline after review, plus the complete formal-suite gate.

Passing TLC cannot close numeric quality, geometry validity, memory magnitude,
latency, renderer semantics, or filesystem crash consistency.  Those remain
independent release evidence even when their control ownership is modeled.
