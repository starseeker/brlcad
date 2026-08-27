# libBObol active debt

Last reviewed: 2026-08-26

This is the authoritative remaining-work list for the Obol drawing stack.
`obol_production_readiness.md` defines the release matrix; `qged_editing.md`
defines the editing workstream; `libbobol_engineering_lessons.md` preserves
resolved failures.  Do not add chronological test logs or alternate designs
here.

## Current baseline

- The semantic drawing API, compact occurrence model, shared-scene/view-local
  split, progressive PoP cache, retained renderer, selection, polygons,
  faceplate, navigation gizmo, and edit-session substrate are implemented.
- **Current architecture gate 2026-08-25:**
  `libbobol_progressive_pipeline_contract.md` is now the canonical ownership
  and complexity contract.  Two duplicate authorities have now been removed:
  the imperative coordinator phase machine (including its event fan-out and
  public diagnostic mirror) was deleted in favor of the derived convergence
  snapshot, and publication batching plus frame acknowledgement now share one
  revision-bound `BObolLodPresentationTransaction`.  The latter atomically
  owns applied-result age/count, frame-request state, reason set, required
  render serial, and retirement.  A failed source also closes its missing-
  inventory edge, so an error cannot remain permanently in discovery without
  another explicit work witness.

  Raw helper-type count is no longer used as the architecture metric: pure
  numeric algorithms and strong value types are desirable, while one mutable
  boolean can still create a competing authority.  Capacity, pressure,
  headroom, point-proxy, and structural facts now share one admission-evidence
  snapshot and pure transition result.  Provider terminality, result age,
  inventory-delta coalescing, and resident-growth obligation now share one
  availability ledger and stateless scheduler; three superseded policy owners
  were deleted.  The admission functional-core boundary is now enforced in
  C++: the controller sees const evidence components, all transitions return a
  complete successor snapshot, and the allocation-free bounded execution
  cursor is physically separate from completed-frame evidence.  Plans and
  cursors now carry the exact inventory, availability, view, policy, and
  capacity revision stamp; stale certified commits are rejected.  Roughly
  forty manual cursor resets were replaced by three typed coordinator revision
  domains plus the existing sole view and policy revision owners.  This closes
  the admission revision gate.  A new
  mutable policy is forbidden unless its unique owner,
  revision token, progress witness, terminal condition, and deleted predecessor
  are identified in the same change.  Preserve measured numeric algorithms and
  compact data structures.  Do not add workload-specific control modes for
  Lucy, multi-Lucy, or 50k/150k.
- **Rechecked 2026-08-26:** after rebuilding every consumer against the current
  Obol headers, the current binary passes all 28 tests carrying the
  `bobol_headless` label in 6.13 seconds real / 20.31 process-seconds,
  including the LoD service,
  cache, compact ownership, retained allocation, edit-manipulator, and GED
  draw-sync contracts.  The draw-sync test now exercises a real independent
  endpoint and proves its local source owner remains isolated across
  draw/erase/redraw/zap.  This is the current headless implementation gate;
  it does not replace the separate graphical release matrix or claim its
  historical aggregate timing/coverage.
- **Rechecked 2026-08-24:** the focused drawing baseline gate passes
  `ged_test_drawing_quad`, `test_qtcad_obol_controller`, and
  `test_qtcad_obol_draw_sync` (44.74 seconds total).  Independent draw
  transactions now bypass the shared scene entirely, preventing a source
  drawn in one independent view from leaking into shared siblings.  The quad
  scene-scope test explicitly suppresses autoview; camera fitting is covered
  separately so asynchronous framing cannot make the scope assertion timing
  dependent.
- Focused TLC models for payload publication, host work, visual arbitration,
  deferred autoview, and Lucy/large-scene convergence pass.  The convergence
  model derives the user-visible modes `interactive`, `discovering`,
  `settling`, `constrained`, and `stable` from controller ownership rather
  than trusting a second mutable HUD state.  On 2026-08-24,
  TLC rechecked the complete model suite without an invariant or liveness
  violation: host work (2,007 generated / 714 distinct), deferred autoview
  (4,833 / 896), arbitration (102 / 72), admission (52 / 40), cold
  presentation (7 / 7), convergence (562 / 355), and spatial production
  (35 / 28).  The convergence model establishes that physical quality debt is
  not itself a quiet residency-work obligation.
  The high-level pipeline model was strengthened on 2026-08-25 to use
  independent inventory, availability, view, policy, and capacity revisions.
  It now proves that a bounded plan interrupted by a superseding availability,
  view, or policy revision can only abort and cannot advance or commit, while
  an unchanged tuple cannot reopen planning.  After adding the retained-
  quality floor across camera epochs, TLC explored 2,986,118 generated /
  1,427,962 distinct states to depth 42 without an invariant, deadlock, or
  liveness failure.
- Exact-current OSMesa shaded interaction lifecycles now pass for true-cold and
  warm Lucy and for warm 50k/150k after the structural-frontier, release-frame,
  deadline-recovery, and coalesced-producer demand-ownership corrections
  described below.  Current warm evidence is the consecutive Lucy pair
  `/tmp/qged-lucy-superseded-fix-{a,b}-20260826` and the 50k/150k pair
  `/tmp/qged-{50k,150k}-superseded-fix-20260826`.  All four finish terminal,
  ready, box-free, and without occurrence failures.  The 150k endpoint is a
  software-performance constraint result, not visual-significance evidence;
  remaining cold, wire, and real-vehicle rows are separate release obligations.
- The coalesced-producer demand contract has a focused formal gate.  On
  2026-08-26 `ObolActiveProducerDemand.tla` passed temporal model checking in
  678 generated / 340 distinct states to depth 12.  Keep this model paired
  with the service/provider regression: the stable asset task, latest demand,
  superseded-result classification, and exact-demand failure lifetime are one
  ownership boundary.  `BOBOL_LOD_PROVIDER_SUPERSEDED` is deliberately distinct
  from genuinely stale cache/source metadata and cannot create an occurrence
  failure.
- **Implemented and qualified on System GL 2026-08-25:** warm spatial PoP
  hierarchies now publish immutable page-local renderer layers instead of
  rebuilding one aggregate mesh.  Lucy's close-view page changes prepare in
  roughly 60--220 ms off the presentation thread, while the former aggregate
  pack alone cost 1.2--1.6 s.  The exact interaction matrix passes: cut 22
  advances to cut 26 during active zoom and cut 36 at rest, zoom-out compacts
  resident memory, the returned view consumes its quality-history proof, and
  no structural boxes or internal page proxies appear.  The obsolete
  append/tombstone/dense aggregate machinery was removed.  Focused cache,
  service, update-action, and live-publication TLC gates pass.
- **Resolved 2026-08-25:** the conservative direct-terminal gate was present
  but never enabled for Generic Twin.  It compared the completed leaf profile
  with a compact presentation count containing the synthetic overview and
  imposed a smaller unique-asset limit even though that profile counts all
  primitive leaves, not only BoTs.  A completed 4.1 MB / 2,248-occurrence
  profile is now admitted under the common 4,096-item bound; per-occurrence
  projection, render cost, worker memory, and resident memory remain the
  authorities.  Cold and warm shaded replays on System GL and OSMesa each end
  with 708 direct full-detail payloads, one progressive payload, all 709 BoT
  occurrences covered, and zero boxes.  qged reports full-detail payloads
  separately and the matrix requires them to dominate this safe scene.
- **Resolved focused OSMesa quality issue 2026-08-25:** provider satisfaction
  and renderer presentation are now independent quality dimensions.  A
  single progressive asset may have its requested suffix resident while a
  motion/quiet-cadence ceiling still hides most of it; that ceiling gap now
  arms the separately bounded, interruptible static-quality pass instead of
  falsely terminating convergence.  The exact-current warm Lucy shaded
  lifecycle restores cut 16 after pose-only rotation, submits 2.55M faces at
  normalized prominent error 2.04, reports zero prominent-floor violations
  and structural boxes, and completes the retained frame in 314--318 ms under
  the 400 ms static-quality bound.  The matrix distinguishes an interrupted
  motion frame from the older completed-frame timing it intentionally
  retains.  Current-binary cold and wire rows remain release work; do not
  generalize this focused warm result to them.
- Current linked-Obol validation also found and fixed an OSMesa retained-wire
  upload overread when a richer GPU prefix survived a coarser request.  The
  150k wire cold/warm replay is now qualified under the bounded diagnostic
  environment.  Discovery now retains one arc
  population and aliases resolved directory names, coverage avoids a second
  whole-file mapping for large databases, and manifest persistence streams
  bounded 8 MiB occurrence chunks before publishing its descriptor, so an
  interrupted write remains a cache miss.  The optional OSMesa flat-wire atlas
  also declines a failed transient allocation and falls back to retained
  rendering rather than throwing through Coin.  Re-run the cold/warm pair to
  establish the resulting end-to-end memory bound and visual result.  The
  compact-presentation path now publishes retained instance records in bounded
  batches rather than creating a second whole-scene update vector; this fixes
  a traced OSMesa repaint `std::bad_alloc` at 150k scale.  Compact discovery
  reserves its lightweight renderer-update population before the first merge,
  avoiding the former 131k-to-262k retained-instance reallocation, while the
  heavyweight occurrence records remain segmented to preserve bounded partial
  progress.  Under the bounded
  8 GiB diagnostic limit, the requested 150k compact population now stops
  cleanly before records are drained, retaining the useful overview/partial
  scene rather than aborting qged.  That limit is below the retained discovery
  footprint of the 3.4 GiB fixture, not another vector-doubling regression.
  Repeat the exact-current cold/warm pairs on adequately provisioned
  production hardware before treating all 150k OSMesa modes as cleared.

  A fresh short cold replay under a 16 GiB address-space cap now discovered
  all 150001 leaves in 20.6 seconds, admitted 664 mesh payloads, and reported
  neither terminal provider failures nor structural boxes.  The constrained
  cache fallback is 512 MiB rather than 4 GiB so it can coexist with a
  multi-gigabyte database under that cap.  This is a regression checkpoint,
  not the required full cold/warm qualification.  The full wire matrix then
  passed cold in 44 seconds and certified warm in 81 seconds, with 150001
  discovered leaves, zero structural boxes or terminal failures, no pending
  work, and an explicit software-performance limit in both terminal frames.
  The corresponding exact-current shaded matrix also passed cold in 44
  seconds and certified warm in 61 seconds with the same coverage, box,
  failure, pending-work, and camera-contract checks.

This is not release clearance.  **Rechecked 2026-08-24:** Generic Twin's
current System-GL and OSMesa shaded/wire cold-warm matrices pass with
`view_ready=true`, zero structural boxes, zero terminal failures, and their
camera contracts in every terminal report.  The exact-current artifacts are
`.build/qged-gui-matrix/20260824-generic-system` (16--21 seconds per
lifecycle) and `.build/qged-gui-matrix/20260824-generic-osmesa`.  The copied
Obol library now has a relative `$ORIGIN` runtime path, so the harness also
verifies that qged uses only `.build/lib` rather than an external renderer
dependency.  The complete System-GL model matrix and visual-importance
quality on real large models remain insufficiently qualified.

## P0: editing runtime contract

- **Resolved CI-registration and coverage gap 2026-08-24:** the existing
  librt primitive edit runtime executables were build-only and consequently
  absent from CTest.  They are now the `edit_runtime` label alongside the
  descriptor contract.  A subsequent audit found that REC and halfspace had
  descriptors but no dedicated runtime executable; REC was also absent from
  the generated internal command-ID header.  The two focused tests now cover
  readback, REC constraint preservation, and rejected REC-radius non-mutation;
  the header is generated from the REC source too.  All 35 runtime tests pass
  in 0.35 seconds.  This preserves continuous coverage of the real primitive
  handlers while the broader advertised-operation rejection and readback audit
  remains open in `qged_editing.md`.

## P0: control-state reduction and maintainability

- **Audit finding 2026-08-25:** the high-level single-pipeline contract and
  exact revision-stamp boundary are sound, but the production C++ shape still
  exposes too much historical control state.  At the start of this audit the
  private coordinator was 953 lines with 40 top-level boolean fields.  Its
  5,930-line policy header held 26 policy/evidence/helper types and another 41
  persistent booleans, while
  controller, progressive-provider, and renderer units contain hundreds of
  direct transition/mutation sites.  These numbers are not themselves proof
  of a defect, but the observed clusters represent overlapping submission,
  retained-allocation, static-trial, point-proxy, interaction-handoff, and
  recovery lifecycles.  Exact revision stamps prevent stale plan commits; they
  do not make all combinations of those latches valid or understandable.
- Treat this as a production P0 gate before further convergence-policy tuning.
  Do not add another Lucy, 150k, cold/warm, renderer, or cache-state branch.
  Those are qualification profiles and numeric planner inputs under the same
  control graph.
- Implement the canonical C++ shape defined in
  `libbobol_progressive_pipeline_contract.md`: one immutable five-domain
  evidence snapshot, at most one revision-stamped bounded plan execution, at
  most one presentation transaction, and a finite typed work ledger.  Outcome
  and HUD state remain derived.
- Migrate behavior-preservingly in deletion-producing slices:

  1. Add the diagnostic-only transition record and offline refinement checker.
     Map every existing production transition to the finite abstract event
     alphabet; an unmapped transition is debt, not a new event by default.
  2. Replace boolean-plus-companion clusters with revision-bound certificate or
     sum types, beginning with static-quality trial and presentation handoff,
     then point-proxy trial, interaction session, retained allocation, and
     submission pass.  Each slice deletes its former fields and reset code.
  3. Introduce one pure reducer which accepts a typed event and returns the
     complete successor control value plus bounded effects.  Move production
     callers behind that boundary; they may no longer write evidence or cycle
     state directly.
  4. Make the controller an effect executor, split evidence/planner/control
     ownership into private units, and delete superseded policy entry points.
     Preserve dense occurrence storage, immutable/moved mesh ownership,
     bounded cursors, and all measured hot-path algorithms.
  5. Add randomized reducer traces and the independent numeric planner oracle
     to CI.  Re-run TLC for ownership/liveness changes, then the complete
     cold/warm System GL and OSMesa profile matrix.  No slice lands on the
     strength of TLA+ or headless tests alone.

- **First deletion-producing slice 2026-08-25:** presentation handoff now uses
  one three-state value (`inactive`, `presentation required`, or `allocation
  required`) in place of three independently mutable booleans.  Two additional
  write-only handoff latches and the always-false
  `preservePresentationLimits` result path were deleted, along with repeated
  reset code and its stale controller commentary.  The existing focused
  handoff lifecycle corpus passes and the value remains trivially copyable.
  This validates the sum-type migration technique; it does not close the
  broader control-state reduction.
- **Second deletion-producing slice 2026-08-25:** static-quality trial state is
  now one allocation-free five-state lifecycle (`idle`, `probing`, `accepted`,
  `reconciling`, or `rejected`) rather than independent active/rejected
  booleans plus a free replay sample.  The type owns its replay sample, resets
  it on every state edge, and refuses to reopen a terminal constraint without
  an explicit reset.  The former `retainAcceptedPhase` helper and its no-op
  controller assignment were deleted.  Focused tests exercise every legal
  state and the rejected-trial non-reopening invariant.  Acceptance and
  rejection are deliberately distinct terminal evidence: a completed bounded
  fractional frame is a retained quality floor, not permission to coarsen as
  though the candidate had missed its capacity bound.
- **Third deletion-producing slice 2026-08-25:** the interaction-start scale,
  scale-valid flag, prior-ready flag, stable budget, and budget-valid flag are
  one allocation-free certificate.  Capture and reset are atomic, so renderer
  or service replacement cannot retain a budget detached from its camera
  evidence.  Focused tests cover empty, inexact-scale, exact-scale, changed-
  scale, and reset states.  This is the evidence certificate for an
  interaction epoch; the remaining interactive/gesture/release phase latches
  still need conversion to one sum type.
- **Fourth deletion-producing slice 2026-08-25:** inventory, availability,
  view, policy, and capacity revision changes now pass through one pure
  `BObolLodRevisionContract`; the public policy-revision synchronization path
  also retires the current cursor.  The same contract classifies an admission
  plan as administrative, exact-current, or stale before commit.  Focused
  executable tests prove that a semantic edge changes exactly one domain and
  that stale certified plans cannot commit.  This is the first production
  refinement-map component, not yet the complete event/obligation trace.
- **Fifth certificate-completeness slice 2026-08-25:** the reversible
  presentation proof now retains fractional next-cut population together with
  integer ceiling, pixel target, point aggregation, population identity, view
  epoch, and measured cost.  The old certificate silently rounded a completed
  `cut N + fraction` frame down to `cut N` when interaction or exact-view
  history restored it.  Focused tests now prove complete-vector pose, quiet,
  and history restoration.  The high-level TLC model was strengthened at the
  same time: a view/input revision cannot replace a retained mesh with a
  coarser stable commit unless a future typed capacity constraint explicitly
  authorizes that regression.  TLC passed 2,986,118 generated / 1,427,962
  distinct states to depth 42.
- **Sixth contract slice 2026-08-25:** renderer-side preparation is now a
  separately modeled finite obligation rather than an unqualified activity
  serial.  `ObolPresentationPreparation.tla` requires one exact target,
  bounded transient reservation, finite total/remaining work, strict
  remaining-work decrease on every retry, and terminal complete/constrained/
  failed evidence when an attempt cannot advance.  TLC explored 1,143
  generated / 526 distinct states to depth 12 without an invariant or
  liveness failure.  The model's first check found an omitted fair admission
  edge, demonstrating the intended guardrail.
- **Seventh deletion-producing slice 2026-08-25:** the production mapping for
  renderer preparation is complete.  Obol's classifier, flat shaded atlas,
  and retained-indirect paths now publish exact-target preparation snapshots
  with an immutable finite total, monotone completed work, retained scratch
  reservation, and complete/constrained/failed terminal state.  libBObol
  aggregates those snapshots at frame boundaries and permits an unchanged
  retry only for a new target, strict same-target progress, or completion.
  The raw `renderPreparationSerial()` API was removed from libBObol's view
  state and controller; Obol retains the counter only as diagnostic telemetry.
  A pre-scan planning obligation is published before consulting an already
  expired host deadline, so an assembly reached late in traversal gets one
  bounded chance to start rather than being falsely classified as capacity
  constrained.  The focused coordinator, update-action, Qt controller, Obol
  CAD, and all 28 `bobol_headless` tests pass.
- **Eighth deletion-producing slice 2026-08-25:** interaction phase, gesture
  kind, last-input time, and settle ownership are one allocation-free
  `BObolLodInteractionSession`.  Callers submit typed begin/input/release
  events and query the session rather than coordinating independent latches.
  `ObolInteractionSession.tla` explored 98 generated / 44 distinct states to
  depth 10, and the focused C++ matrix covers the legal session transitions.
- **Ninth identity-hardening slice 2026-08-25:** retained allocation now has
  one canonical semantic input key shared by reuse certificates and in-flight
  transaction matching.  The disabled maximum-protected-budget value is
  normalized out of the key; changing that dormant raw input no longer
  advances the plan serial.  The focused update-action test exercises the
  inactive-input case.
- **Tenth ownership slice 2026-08-25:** a render deadline now selects one typed
  successor in priority order: strict finite-transaction retry, continuation
  of an existing population cursor, or presentation-capacity recovery.  The
  old path could reset a live 150k cursor to zero, advance the capacity
  revision, and request another frame in the same event.  The four-case C++
  matrix and `ObolDeadlineOwnership.tla` guard the replacement; TLC passed 43
  generated / 34 distinct states to depth 5.  The first full 50k replay then
  exposed an incomplete production refinement mapping: population continuation
  preserved a replay latch left by the former render owner, so the progressive
  pump returned before reaching the live cursor.  The effect now retires that
  latch, and a focused controller test reproduces the owner transfer.  The
  same full warm OSMesa script which failed after 194 seconds now passes in 42
  seconds, terminal and box-free, including motion, zoom, selection, subpath
  erase, and redraw.
- **Eleventh deletion-producing slice 2026-08-25:** quiet restoration now has
  one allocation-free `BObolLodQuietSuccessorReducer`.  It selects the complete
  target and typed handoff from revision-matched exact-view, prior-stable,
  proven-scale, or stable-demand evidence using fixed precedence.  A restored
  prior pose must first produce one current-pose presentation proof; transient
  motion limits otherwise request at most one bounded allocator handoff.
  Their arrival order cannot select the semantic quiet target.  The
  controller's later exact-history overwrite and separate one-pixel fallback
  were deleted, as was the older release-time pose restore which could publish
  capacity evidence before the quiet transaction owned the target.  Release
  now retains the responsive motion image through debounce; only the quiet
  reducer restores the prior pose.  A focused C++
  permutation test holds semantic evidence fixed while varying every transient
  renderer control, and `ObolQuietSuccessor.tla` passed 1,242 generated / 924
  distinct states to depth 5 with schedule independence, current-pose proof,
  and liveness intact.
  Nested erase/redraw was hardened in the same qualification pass: changing a
  retained occurrence-visibility frontier no longer advances mesh-inventory
  demand or evicts its payload.  Headless Qt tests prove inventory and payload
  identity remain unchanged while presentation revisions advance, and the
  graphical Lucy redraw now returns directly to the retained mesh rather than
  exposing an empty intermediate frame.
- **Twelfth control slice 2026-08-25:** the exact-current 150k OSMesa trace
  proved the
  source census itself is no longer the long pole: 2,048-entry windows took
  about 5.7--7 ms and a continuously pumped full pass took less than a second.
  The retained pixel-demand endpoint now uses a revision-bound capacity-search
  certificate with ordered scene-wide candidates, three samples per candidate,
  and at most one steady-to-static goal transition.  A frozen search must
  consume a sample, narrow its safe/unsafe bracket, or publish terminal typed
  evidence.  `ObolCapacitySearch.tla` checks every modeled steady/static
  capacity combination (1,792 generated/distinct states, depth 31), and the
  focused C++ matrix covers safe, unsafe, stale, and unmeasurable samples.
  Object count and workload identity do not select this path.  The migration
  is not complete: ordinary coverage calibration still has legacy probe/count/
  rescan fields.  Move those consumers behind the same certificate and delete
  their old authority rather than layering another search over it.  Current
  150k System-GL/OSMesa qualification remains a production P0.
- **Thirteenth deletion-producing slice 2026-08-25:** quiet point calibration
  and retained triangle recovery are now one allocation-free three-state
  `BObolLodPointQualityPhase` instead of independently mutable pending flags.
  Recovery has priority because it needs source admission; calibration pauses
  that producer.  The former combination could advertise both phases and leave
  a 50k scene pending with no enabled task, submission, frame, or render.
  A second replay exposed the production refinement seam: a completed recovery
  pass changed its cut, generic cleanup cleared the persistent changed latch,
  and no presentation frame was scheduled.  The sum type now consumes the
  immutable completed-pass result to require that frame.  Focused C++ tests and
  `ObolPointQualityOwnership.tla` cover the ownership transfer; TLC explored 6
  generated / 5 distinct states to depth 4 without an invariant, deadlock, or
  liveness error.  Named retained recovery, handoff, reallocation, and rescan
  obligations now also enter the derived convergence outcome through one
  `controlPending` input; they can no longer coexist with a reported terminal
  view merely because no source task or frame is currently active.
- **Fourteenth liveness slice 2026-08-25:** rejected static-quality
  reconciliation now has one terminal transition and one allocation identity.
  The old path left a renderer ceiling after occurrence-local handoff, then
  changed an irrelevant throughput-derived budget in the retained-allocation
  key, so the same fixed presentation reopened indefinitely.  A second trace
  found two deadline owners: the terminal rejected trial required the 400 ms
  static bound, but older structural-repair state replaced it with a 200 ms
  calibration bound and restarted capacity recovery.  Reconciliation now
  consumes its one-shot state, normalizes the fixed presentation budget, and
  terminal/static ownership wins deadline selection.  The exact warm 50k
  OSMesa lifecycle passes in 93 seconds with 7,267 mesh payloads, 42,733
  subpixel occurrences, zero boxes, zero foreground obligations, and a held
  terminal result after motion.
- **Fifteenth control-boundary slice 2026-08-25:** production convergence and
  diagnostics now consume one allocation-free refinement map from the
  remaining concrete latches to nine finite obligations and one selected
  owner.  The former duplicated calibration/control Boolean unions were
  deleted.  The map rejects ownerless work, an owner which does not own a
  present obligation, terminal foreground work, and invalid readiness.  Its
  complete 20-input Boolean domain (1,048,576 combinations) is executable in
  the focused C++ test in 0.14 seconds; qged requires the obligation, owner,
  and zero-violation fields on every graphical sample.  Exact warm OSMesa
  qualification passes at both 50k (93 seconds) and 150k (68 seconds); the
  latter holds a zero-obligation terminal state with 3,560 mesh payloads and
  no boxes after its complete interaction/selection lifecycle.  This is the
  first live refinement boundary, not yet the final event/effect reducer:
  subsequent slices must move effect writers behind it and delete their old
  entry points.  `ObolControlRefinement.tla` independently checks both error
  values across all 1,048,576 fact combinations.  TLC generated 4,194,302
  states / 2,097,152 distinct states without an invariant, liveness, or
  deadlock error; fair owner service strictly decreases the concrete fact set
  and eventually retires all work.
- **Sixteenth successor-ownership slice 2026-08-25:** submission activity and
  full-rescan debt are one `BObolLodSubmissionPass`, and numeric capacity
  sampling owns its complete unchanged-frame sequence.  Retiring the generic
  presentation barrier for the same framebuffer may not also start an
  occurrence pass.  The old dual effect inserted one O(scene) allocation
  between every capacity sample and made a 150k view visibly alternate between
  refining and balancing.  `ObolSubmissionPass.tla` covers paused discovery
  debt; the capacity-search model covers the exclusive sample sequence.
- **Seventeenth resident-growth slice 2026-08-25:** immutable suffix discovery
  now has one typed availability-ledger phase: idle, drain required, drain
  active, dirty active drain, or reallocation ready.  The separate controller
  drain boolean was deleted.  Completing a drain can only require another
  coalesced drain or enable one retained reallocation; it cannot fall through
  to ordinary capacity calibration.  The prior fallthrough repeatedly ran an
  unchanged 150k suffix pass because allocation was correctly prohibited while
  the resident population was unsettled.  `ObolResidentGrowth.tla` checked 155
  generated / 119 distinct states to depth 19, including growth arriving
  during a drain and eventual quiescence.  The focused ledger test covers the
  same transitions and interruption by a new view/source epoch.
- **Eighteenth point-recovery liveness slice 2026-08-25:** point-to-triangle
  recovery completion is one idempotent controller operation shared by pass,
  completed-frame, and progressive-pump effects.  A 150k OSMesa run reached
  99 percent with no submission, worker, result, or render work, but retained a
  recovery phase because a stronger barrier had consumed its frame callback.
  The already presented recovery now retains a level-triggered pump witness.
  The updated `ObolPointQualityOwnership.tla` checks both direct-frame and
  deferred-pump retirement (8 generated / 6 distinct states, depth 5).
  Exact-current graphical evidence at
  `/tmp/qged-control-terminal-150k-final-20260825` reaches a 100-percent,
  box-free, zero-owner terminal state after the complete interaction/selection
  lifecycle.  Its first in-run validation used the preceding 280 ms stable
  framebuffer as the held-motion duration even though held frames were
  deliberately interrupted against their 100+5 ms deadline.  The specialized
  assertion now uses that typed interruption witness.  It also uses the current
  heterogeneous render-cost budget rather than the removed face-budget field;
  a quiet unchanged image may use the separate 400 ms static-quality allowance.
  Revalidation of the captured report and images passes the complete matrix
  validator (`validation-rerun2.log`).  Lucy System GL remains green after the
  shared completion extraction at
  `/tmp/qged-point-recovery-lucy-final-20260825`.
- **Nineteenth HUD-ownership slice 2026-08-25:** periodic LoD telemetry no
  longer calls the general GED faceplate synchronizer.  The new narrow
  `ged_view_lod_progress_sync` entry point publishes only the retained LoD
  track, fill, and label records; axes, grid, ADC, framebuffer composition,
  and unrelated HUD state remain untouched.  The focused faceplate test
  stages a grid change and proves progress-only publication preserves the
  existing grid node while a later full sync applies the change.  The GED API
  manifest, faceplate drawing test, qtcad faceplate test, and exact qged
  Generic Twin OSMesa shaded cold/warm lifecycle pass; its report and images
  are `/tmp/qged-lod-hud-isolation-20260825`.  Retain this boundary:
  progress text is an observer and may never manufacture general scene
  mutation or a convergence successor.
- **Twentieth trace-refinement slice 2026-08-25:** the public convergence
  snapshot and qged report now carry the complete five-domain revision stamp,
  retained allocation-plan identity, presentation-transaction identity and
  required render serial, and committed-frame serial.  The GUI matrix runs a
  standalone jq refinement checker over every controller sample.  It rejects
  revision/transaction/frame regression, multiple nonzero plans certified by
  one unchanged revision tuple, spontaneous terminal reopening, nonzero local
  invariant masks, observer disagreement, and an A/B/A control cycle with no
  changed revision, plan, transaction, frame, or input witness.  Repeated
  identical states remain legal because bounded worker/cache activity may
  legitimately outlive a sample interval.  A focused checker test proves
  valid held work and explicit command reopening are accepted while synthetic
  plan replacement, terminal reopening, and A/B/A traces fail.  The exact
  Generic Twin OSMesa shaded cold/warm lifecycle passes at
  `/tmp/qged-control-trace-final-20260825`, ending idle and terminal with no
  control owner or obligation; no qged process remains.  This is a sampled,
  observational acceptance boundary and never participates in scheduling.
  Recording every production transition with its abstract event and finite
  remaining-work measure remains part of the canonical reducer migration.
- **Twenty-first availability-successor slice 2026-08-25:** an exact warm 50k
  trace found the remaining pre-drain form of the resident-growth cycle.  A
  completed ordinary pass observed a new resident suffix, requested capacity
  calibration, and restarted the shared submission cursor; the availability
  scheduler could not begin its required drain until that cursor became idle.
  `BObolLodAvailabilityScheduler::completedPassSuccessor` now selects one pure
  successor with drain completion before resident-growth yield and capacity
  calibration last.  The controller applies that decision once and no longer
  lets capacity manufacture a pass while availability owns the population.
  The focused permutation test covers the precedence, and
  `ObolResidentGrowth.tla` now includes the in-flight ordinary cursor and
  rejects capacity restart after a growth edge.  TLC checked 155 generated /
  119 distinct states to depth 19 with no invariant, liveness, or deadlock
  error.  Exact OSMesa shaded 50k qualification passes cold in 47 seconds and
  warm in 56 seconds, including motion, zoom, hierarchy, erase/redraw, and
  selection; both terminal views have zero boxes and no foreground owner.
  The same current binary passes the 150k warm OSMesa lifecycle in 91 seconds;
  its initial and post-motion waits terminate in 48.5 and 19.5 seconds with
  zero boxes, zero foreground obligations, and a final 4,253 mesh / 147,933
  aggregate-point presentation.
  Repeated force-terminal pumping was also made edge-triggered, and structural
  preloading now uses the renderer's full point-proxy threshold rather than a
  second incompatible hysteresis boundary.
- **Applied-allocation floor ownership 2026-08-26:** the warm 50k OSMesa wire
  lifecycle exposed a fully current occurrence plan whose selected cost was
  only 0.14 percent above its certified budget.  The generic capacity fallback
  nevertheless treated that completed result as missing evidence and rebuilt
  it indefinitely.  `claimOverBudgetAllocation` now assigns this state to the
  presentation-reconciliation handoff, whose only successors are local
  small-occurrence aggregation, one bounded static attempt, or a terminal
  constraint.  The focused coordinator regression passes and
  `ObolCapacityPresentationHandoff.tla` checks 178 generated / 138 distinct
  states, including both reducible and irreducible floors.  The exact-current
  untraced GUI lifecycle now passes in 74 seconds with no boxes, owner, or
  obligation, 3,210 mesh payloads, 1.02M faces, 254 capacity revisions, and 53
  plans.  The prior untraced run took 133 seconds; the diagnostic replay still
  cycled after 220 seconds.
- **Timer-edge handoff evidence 2026-08-26:** the warm Lucy OSMesa wire zoom
  completed its corrected renderer frame at roughly the 100 ms endpoint.  A
  frame just over that edge correctly supplied no new safe extrapolation, but
  `noteFramePresented` also discarded the stricter 1.3M-work recovery limit
  recorded by the interrupted richer frame.  The handoff then accepted a
  3.38M-work plan, removed its ceiling, and retried the same known deadline
  miss indefinitely.  Presentation completion now preserves the prior limit
  when no safer replacement sample exists.  The coordinator regression covers
  that boundary, and the exact untraced lifecycle passes in 68 seconds through
  smooth zoom, rotation, hierarchy selection, and subpath erase/redraw.  Its
  final 2.10M-face view is terminal with no boxes, owner, or obligation and a
  0.632-pixel error certificate.
- **Availability-bounded allocation and truthful HUD 2026-08-26:** a warm
  System-GL Lucy rotation reached a current memory denial at cut 26 while its
  retained allocator repeatedly certified unavailable cut 27.  The provider
  correctly suppressed the suffix, but the presentation handoff waited for an
  impossible assignment at 75 percent.  Retained-allocation identity now
  includes the resident-admission revision.  A current denial restricts only
  the allocation endpoint to the active resident cut; the certificate retains
  the full pixel-demand cost so reclaimed capacity can reopen richer work.
  Applying that same-or-coarser resident allocation preserves the denial
  witness and performs no provider I/O.  The focused action test covers this
  refinement, and `ObolCapacityPresentationHandoff.tla` now requires assigned
  cuts to become drawable before application (178 generated / 138 distinct
  states, depth 9).

  The same replay exposed a separate observation race: the controller had
  entered interactive state, but the retained label still said `View ready`
  until the next paint.  Beginning and ending a pointer interaction now
  publishes the cheap progress-record transition synchronously while leaving
  rendering queued and coalescible.  The GUI report records the actual label
  text and rejects `View ready` unless convergence has state, is terminal and
  view-ready, has no background work, and is at 100 percent.  Targeted warm
  shaded Lucy lifecycles pass this contract on System GL in 16 seconds and
  OSMesa in 74 seconds at
  `/tmp/qged-lucy-system-hud-interaction-fix-20260826` and
  `/tmp/qged-lucy-osmesa-hud-interaction-fix-20260826`.

  A later Lucy checkpoint exposed the corresponding non-pointer edge: a
  completed terminal frame could schedule background compaction on the host
  wake while the retained label still described the preceding state.  The Qt
  frame-request callback now synchronizes the progress-only faceplate
  immediately after its bounded pump and before deciding whether a geometry
  frame is required.  `/tmp/qged-hud-wake-lucy-20260826` passes the full warm
  OSMesa lifecycle in 76 seconds with zero ready-label/state mismatches.  The
  50k and 150k OSMesa replays below have the same zero-mismatch result.
- **Observability debt made explicit 2026-08-25:** the GUI matrix contained
  three coordinator-invariant assertions using jq's `// 0` fallback, but qged
  never emitted those fields.  They therefore passed by construction and were
  removed.  Do not restore placeholder telemetry.  The sampled refinement
  checker described above now requires real fields with `has(...)`; it does
  not substitute for the remaining per-transition abstract-event record.
- **Foreground/background qualification correction 2026-08-25:** a warm 50k
  OSMesa view reached terminal constrained presentation with zero boxes and no
  source, admission, or frame owner while thousands of finite resident-prefix
  trims continued in the service.  The interaction matrix formerly called its
  strict cache-idle wait at that point and reported a 180-second visual
  convergence failure even though the view was ready.  50k/150k interaction
  checkpoints now wait for `viewReady`; background compaction retains its own
  explicit completion/memory qualification.  Do not make cache persistence or
  resident normalization a prerequisite for accepting input.  The corrected
  warm shaded 50k OSMesa lifecycle passes in 62 seconds, including motion and
  zoom while compaction remains active, with terminal constrained coverage,
  7,252 mesh payloads, zero structural boxes, and no foreground owner.  Its
  retained evidence is `/tmp/qged-contract-50k-visible-ready-20260825` for the
  duration of this qualification session.
- **Renderer reservation conformance closed 2026-08-26:** resumable classifier
  scratch is now dense, capacity-accounted storage; only the sparse published
  change sets retain hash containers.  Per-part frame-plan bookkeeping shares
  the same dense index domain instead of allocating several independent hash
  tables.  Obol reports every retained scratch capacity through
  `reservedBytes`, and its CAD regression exercises a 131,072-occurrence
  preparation.  This closes the known unaccounted classifier reservation; any
  future resumable container must extend the same certificate before it may be
  admitted.
- **Closed common-contract counterexample 2026-08-25:** typed preparation alone
  left quiet camera restoration split between the ordinary quiet-capacity
  pass, exact-view history, and a later static-quality trial.  Diagnostic
  timing could consequently change the terminal Lucy cut.  The eleventh slice
  above moved selection of the first quiet target behind one reducer and
  deleted the two downstream overwrite paths.  Future quality tuning must
  consume that target through the common capacity-search contract; it may not
  recreate a pose, Lucy, OSMesa, trace-enabled, or object-count exception.
- **Rejected common-restore experiment 2026-08-25:** immediately entering a
  static-quality trial for every certified restored presentation improved the
  terminal Lucy pose cut from 14 to 20, but it was not a common solution.  A
  50k warm OSMesa replay failed to reach its first stable scope within 180
  seconds and a 150k cold replay failed to settle after rotation within the
  same bound.  The production change and its speculative helper/tests were
  reverted; all 28 `bobol_headless` tests pass on the restored baseline.
  Preserve this as a design counterexample: target restoration and permission
  to create new static-quality work are separate reducer events.  Do not add a
  cardinality branch to distinguish them.
- **Closed prominent-floor counterexample 2026-08-25:** an earlier warm shaded
  OSMesa Lucy replay terminally constrained at cut 16 with one protected-floor
  violation despite the 400 ms static allowance.  The current retained-floor,
  complete-certificate, adaptive spatial range, and quiet-successor paths
  reach the close, zoom-out, and exact-return checkpoints with zero protected-
  floor violations and no structural boxes.  The zoom-out may intentionally
  use the event-driven static-quality allowance; requiring the ordinary redraw
  cadence for that unchanging image was an obsolete test assertion, not a
  reason to discard affordable fidelity.
- **Qualification-contract correction 2026-08-25:** the qged matrix now
  distinguishes a useful cold structural preview, an adopted progressive CAD
  mesh, current-view terminality, and optional background cache persistence.
  A cold Lucy smooth-zoom test formerly began against the first structural
  preview and later required global background work to be idle even though
  the foreground view was ready.  The corrected System-GL shaded cold/warm
  pair passes in 112/21 seconds with cut growth, resident compaction, exact-
  view history recall, zero boxes, and no prominent-floor violations.  Keep
  these four facts separate in clients and the HUD.
- **Closed recovery/capacity ownership counterexample 2026-08-26:** an exact
  warm 50k OSMesa interaction replay could remain at 95--99 percent for more
  than 180 seconds with no worker or result work.  A point-to-triangle
  recovery allocation had reached its finite certified population, but its
  intentionally unsatisfied richer pixel demand was mistaken for incomplete
  admission.  Ordinary capacity calibration then restarted the same compact
  occurrence scan; failed runs advanced 950--3,017 capacity revisions without
  a terminal witness.  The availability successor now gives point recovery
  exclusive capacity ownership after resident growth, and recovery accepts
  residual quality debt only when the current allocation certificate covers
  the current population.  The focused coordinator test covers both sides of
  that proof.  `ObolPointQualityOwnership.tla` checks exclusive capacity/point
  ownership, changed-frame completion, deferred completion, and certified
  no-op completion in 13 generated / 8 distinct states.  Three untraced 50k
  OSMesa replays now terminate at 100 percent in 28--53 seconds with all
  active payloads satisfied and 112--183 capacity revisions; the same binary
  terminates the 150k OSMesa rotation replay at 100 percent with no owner,
  obligation, boxes, or background work.
- **Point-frame ownership completion 2026-08-26:** a later warm 50k replay
  exposed three adjacent ambiguities rather than a new workload regime.  A
  rejected one-pixel static trial retained the point-calibration owner and
  retried its rejected framebuffer; a point request arriving during an
  already-active finite capacity search allowed that search to reproduce
  itself before the point frame; and the one phase named both adaptive
  point-cut calibration and confirmation of an already chosen handoff cut.
  Static rejection now restores its recorded retained threshold and releases
  both trial owners atomically.  Completed-pass scheduling gives a waiting
  point frame precedence over starting another capacity search.  The point
  quality sum type now distinguishes `ADAPTIVE_CALIBRATION`,
  `HANDOFF_CONFIRMATION`, and `TRIANGLE_RECOVERY`, so confirmation cannot
  immediately undo the cut it was meant to prove.

  The same maintainability pass replaced retained allocation's pending,
  preserve-budget, and reconciliation-budget field trio with one finite
  request value; invalid companion-field combinations are no longer
  representable.  Focused tests cover every alternative and successor edge.
  `ObolPointQualityOwnership.tla` now covers static rejection, capacity/point
  overlap, handoff replacement, and recovery in 29 generated / 17 distinct
  states to depth 5 with no error; TLC found and forced correction of a missing
  capacity-drain guard before the production change was accepted.  The final
  headless label passes all 28 tests.  Exact warm OSMesa Lucy and 150k
  interaction lifecycles pass.  The 50k lifecycle no longer cycles: it either
  reaches a terminal ownerless view or shows finite advancing plan work.

  **Renderer/deadline follow-up resolved 2026-08-26:** `calibrationPending`
  still grouped adaptive point census with confirmation of an already chosen
  handoff cut.  The latter therefore inherited a census-only 400 ms allowance
  and completed a 342.87 ms post-release OSMesa frame.  The broad predicate is
  now named `presentationPending`; only `adaptiveCalibrationPending`
  authorizes that extension.  A handoff may still use the independently
  justified static-quality deadline when its allocation certificate requires
  it.  OSMesa's
  synchronous flat-shaded executor also polls cancellation every 8K rather
  than 64K triangles; perf showed software lighting and rasterization, not
  controller work, dominate this endpoint.  At
  `/tmp/qged-handoff-deadline-50k-20260826` the over-budget post-release
  candidate is interrupted at 122.39 ms, the previous complete 254K-triangle
  framebuffer remains visible, and the run passes in 52 seconds before
  terminating at 1.29M faces with no box or HUD mismatch.  Do not restore the
  census extension to handoff confirmation merely because it is a point
  presentation, or enlarge software draw chunks beyond the host's
  cancellation granularity.
- **Interactive cost-domain recovery closed 2026-08-26:** the next 50k trace
  showed that the renderer honored its 100 ms motion deadline, but the
  controller reduced its ceiling by only one PoP ordinal after every miss.
  Several over-budget frames could therefore expire without ever presenting
  the changed camera pose.  Interactive and quiet deadline misses now use the
  same measured attempted-cost/deadline ratio to choose the next ceiling in
  render-cost space.  Only quiet recovery may retain the prior-pose
  deadline-safe floor; changed-pose motion has no such proof.  The current
  warm OSMesa replay reaches ceiling zero directly, presents the held-motion
  pose in approximately 65 ms, and then reaches its zero-box 1.52M-face static
  endpoint in 78 seconds.  This is a common numeric recovery rule, not a 50k
  mode or object-count threshold.
- **Structural-frontier liveness and release ownership closed 2026-08-26:**
  an exact warm 150k OSMesa frame reached 79 percent with 5,160 structural
  boxes while every provider, task, result, submission, and control owner was
  idle.  The point-quality model did not include this renderer-owned frontier,
  so a completed capacity sample could retire without choosing its successor.
  The derived control mapping now treats an unresolved structural frontier as
  planning work.  When its exact point distribution has another finite
  threshold, that frame owns the successor; otherwise ownership transfers to
  bounded structural repair.  No whole-scene source rescan is introduced.
  `ObolStructuralFrontierOwnership.tla` exhaustively checks 30 distinct states
  to depth 14, including the invariant that every ownerless nonempty frontier
  has an enabled successor and the liveness property that finite work reaches
  readiness.

  The same full interaction replay found a second, independent release edge.
  `endLodInteraction()` documented retention of the responsive motion
  presentation through debounce but immediately disabled renderer camera-frame
  reuse.  That forced an exact 150k reclassification and completed a 410 ms
  release frame.  Release now leaves reuse owned by the typed interaction
  session; the quiet-successor transition disables it and begins exact
  classification after debounce.  The final warm OSMesa matrices pass for
  Lucy, 50k, and 150k.  The 150k run retains its 142 ms completed motion frame,
  records a bounded 114 ms interruption without advancing the completion
  serial, then reaches a box-free ownerless terminal view in 95 seconds.
  The GUI validator accepts an interrupted motion witness only when the
  completion serial is unchanged and the interruption met its deadline; a
  newly completed over-budget frame still fails.
- **Submission source ownership reduction 2026-08-26:** the bounded submission
  cursor formerly kept source identity, entry population, validity, and a
  write-only retained-admission marker as independent fields.  Its non-compact
  visibility census likewise kept a source pointer and count which callers had
  to reset in lockstep.  They are now two small owner-thread values: source
  identity is the plan's validity witness, and changing the census source
  atomically starts a new saturating count.  The dead marker and repeated reset
  code are gone; no planner or scheduling branch was added.  The dependency-
  complete 28-test headless gate passes.  Exact-current warm OSMesa lifecycle
  replays also pass at `/tmp/qged-submission-owner-lucy-20260826` and
  `/tmp/qged-submission-owner-50k-20260826`: Lucy terminates ready with 1.95M
  presented faces, 0.632-pixel maximum certified error, and no boxes or
  failures; the 50k compact scene terminates ready at 1.32M faces with no boxes,
  failures, or pending work.  This is a state-representation improvement, not
  another workload policy.
- **Submission delta ownership reduction 2026-08-26:** exact inventory and
  structural-repair deltas formerly stored an independent active latch, target
  source vector, and optional selective-plan vector.  Those fields permitted
  active-without-target, plan-without-target, and stale-plan combinations.
  One typed value now derives activity from its target population and
  represents each target as either a full source scan or a selective dense-
  entry plan.  Widening selective work to a full scan is explicit, and an
  already-full target cannot accidentally be narrowed by a later append.
  Repeated reset/search/erase code and the independent active latch are gone.
  A focused value test covers initial/reset state, duplicate targets,
  selective-to-full widening, non-narrowing, and replacement-plan merging.
  All 28 headless tests pass.  The exact-current warm OSMesa Lucy and 50k
  lifecycles at `/tmp/qged-submission-delta-{lucy,50k}-20260826` also pass:
  both terminate ready with zero structural boxes, occurrence failures, or
  pending work; their final certified errors are 1.265 and 3.000 pixels.
- **Structural-repair transaction ownership 2026-08-26:** the exact fallback
  frontier was represented by a pending latch plus independently reset
  frontier-count and cost-reservation fields.  It is now one value whose
  activity is derived from a nonzero exact frontier.  Beginning or retiring
  the frontier clears its reservation, a reservation cannot attach to an
  inactive frontier, and repeated budget evaluation cannot silently replace
  the first reservation for the transaction.  The stale point-recovery latch
  comment left behind by an earlier extraction was removed as well.  Focused
  tests cover inactive admission, begin, one-shot reservation, reservation
  reset, empty frontier, and retirement; all 28 headless tests pass.  The
  canonical exact-current warm 150k OSMesa regression at
  `/tmp/qged-structural-repair-owner-150k-20260826` passes in 101 seconds and
  terminates ready at 1.74M faces and 1.000-pixel certified error, with zero
  structural boxes, failures, background work, or other pending control work.
- **Optional-value cleanup 2026-08-26:** forced PoP cut and exact camera
  identity each used a value plus an independently mutable validity boolean.
  Forced cut is now `std::optional<int>` and exact camera identity is an
  optional snapshot; their absent states cannot expose stale companion data.
  Existing forced-cut integration coverage and all 28 headless tests pass.
  The exact-current Generic Twin shaded cold/warm matrix at
  `/tmp/qged-optional-policy-signature-generic-20260826` passes on System GL
  in 12/11 seconds and OSMesa in 15/14 seconds.  Every terminal frame is
  ready and box/failure-free, every lifecycle records one exact-view history
  recall, and the paired final images have the same camera and scene content.
- **Spatial residency observation correction 2026-08-25:** a global PoP cut
  does not order memory use for a chunked mesh.  A wide view may need hundreds
  of coarse pages and legitimately retain more bytes than a clipped close view
  using fewer rich pages.  The service now publishes a per-consumer certificate
  for the last complete demand-aware compaction scan: exact demand revision,
  resident-generation currency, and real trim/eviction candidate count.  The
  Lucy zoom test accepts either measured byte reclamation or a current
  zero-candidate certificate under the configured memory limit; it no longer
  forces useful visible pages out to make a global byte counter fall.  The
  current System-GL shaded lifecycle passes in 16 seconds.  One observed run
  reclaimed the 321 MiB close-view peak to 258 MiB on zoom-out with one real
  compaction; foreground readiness remained independent of the explicitly
  owned background scan.

- Acceptance is structural as well as behavioral: the finite event alphabet
  and state variants are enumerable from the control header; no renderer, HUD,
  timer, or client callback owns retry/readiness policy; unchanged evidence
  cannot reopen planning; and every removed latch remains absent.  The full
  graphical matrix must show no regression in first useful image, convergence,
  visual importance, FPS, memory recovery, selection, or editing.

## P0: visual quality and progressive convergence

- **Resolved 2026-08-23:** cache-independent deferred-autoview ownership.  The
  2026-08-23 OSMesa
  Generic Twin shaded cold/warm real-canvas pair reached the same 709 mesh
  payloads with zero boxes, failures, or pending work, but its `ae 90 0`
  stable camera differed materially by cache state.  Cold deferred autoview
  and warm cached completion must target the same view context and preserve
  the same caller orientation before release.  The retained evidence bundle
  is `/tmp/qged-generic-smoke-20260823-1947` until the corrected replay is
  captured; do not weaken the camera-contract comparator to hide it.  The
  immediately repeated pair passed, including when timing diagnostics slowed
  execution, so treat this as a timing-sensitive synchronization race rather
  than a deterministic cache-format discrepancy.
  The retained failing telemetry identifies the concrete software-canvas
  asymmetry: cold reaches the requested `ae 90 0` camera after its first
  progressive synchronization, while warm remains at the initial camera.
  `QgSW::need_update` previously requested only a framebuffer repaint after
  `diff_hashes()`, unlike `QgGL`, so it did not copy a GED command's camera
  into the Obol controller.  The source now requests `BV_REFRESH_VIEW` before
  presenting.  A rebuilt qged passed the exact cold/warm OSMesa replay in
  `/tmp/qged-generic-camera-fix-20260823-rerun`: both runs recorded the same
  `ae 90 0` orientation and the matrix reported zero failures.  Retain the
  old failure bundle as a regression specimen and preserve this replay in the
  release evidence until the full GUI matrix supersedes it.
  The first fix exposed the same ownership edge in a fast warm Lucy wire
  presentation: the first background frame could be requested without a
  semantic canvas refresh.  `QgSW::present_frame` now includes
  `BV_REFRESH_VIEW`, which synchronizes the authoritative GED view at the
  presentation boundary without forcing every worker wake to mutate camera
  state.  The validated warm replay in
  `/tmp/qged-lucy-wire-present-fix-20260823` completed in 61 seconds and its
  stable camera exactly matched the preceding cold wire contract.
- `ObolDeferredAutoview.tla` is the companion bounded ownership model for
  that investigation.  It verifies request replacement, same-view framing
  cancellation, orientation preservation, exact-once completion, teardown
  without cross-view writes, and the producer/wake/presentation liveness
  witness.  TLC passed this extension on 2026-08-23 (4,833 generated / 896
  distinct states, depth 11).  It does not validate Qt event ordering or
  floating-point camera values; retain the real qged cold/warm camera replay
  as the implementation-level acceptance test.
- **Resolved 2026-08-24:** normal orthographic zoom no longer narrows the
  camera depth range.  The direct camera now retains the largest view scale
  for the current scene structure as its depth reference, so a close zoom
  changes only projected scale rather than becoming an implicit section cut.
  `libBObol_lod_update_action` covers this explicitly.  A System GL Hubble
  replay drove the view from 1176.54 mm to 0.0001 mm while retaining focal,
  near, and far distances of 5882.72, 1.1765, and 123537.19 respectively,
  with zero fallback boxes.
- **Resolved 2026-08-24:** a resize replay exposed a second autoview
  violation: provisional coverage bounds could satisfy the deferred request,
  then a later exact bound recentered the same view during fullscreen.  Only
  completed producers and verified exact root proxies may now fulfill the
  request.  The formerly failing Generic Twin System GL resize row passes
  with identical center and scale before fullscreen, during fullscreen, and
  after restoration.  Requalify the complete resize matrix after this change.
- **Resolved 2026-08-24:** switching a settled retained draw from `view lod 1`
  to `view lod 0` cleared its view payloads and exposed structural boxes until
  another realization occurred.  Disabling adaptive selection now preserves
  the current immutable presentation; re-enabling starts a new policy epoch.
  Generic Twin `-m1` resize/policy-toggle replays pass on both System GL and
  OSMesa with no fallback boxes, stable camera, and correct HUD state.
- **Resolved 2026-08-24:** OSMesa hidden-line (`draw -m4`) on Generic Twin
  has a terminal structural-coverage path.  A compact
  view previously charged its startup boxes against the only initial mesh
  allowance and could therefore remain at 673 boxes forever.  Point
  aggregation now cannot preempt source admission before an assembly exists,
  and a cold view may spend one bounded 500k-cost replacement wave.  Once
  point aggregation reaches its maximum threshold, the remaining boxes use a
  one-shot structural repair transaction rather than cycling its calibration
  latch.  The direct OSMesa smoke now reaches 244 CAD payloads and zero
  visible structural boxes.  A one-shot coverage frame may use a 200 ms
  deadline, and the static framebuffer handoff may complete without an
  interactive deadline; both are bounded quiet-state exceptions.

  A broad point-aggregation eligibility check previously vetoed the
  irreducible-frame escape even after every aggregation transition was
  settled, repeatedly aborting the same 101 ms OSMesa frame.  The controller
  now distinguishes an armed aggregation transition from mere eligibility.
  The full OSMesa m4 Generic Twin resize row passes in 21 seconds: its stable
  checkpoint has exact CAD work, 286 payloads, zero structural boxes, a 1.0
  convergence fraction, no stable handoff, and `view_ready=true`.  Retain the
  zero-box and terminal-idle assertions in the full matrix.
- **Resolved 2026-08-24:** a cold Hubble System-GL draw could either omit its
  root overview or retain it permanently.  The overview is deliberately
  outside leaf replacement accounting, but is now a visible retained
  structural primitive in shaded and wire modes.  It remains until the
  retained assembly has committed useful leaf coverage; a complete leaf-box
  frontier is sufficient, since the overview is only an early extent cue, not
  required geometry detail.  Detached/headless sources retire it at complete
  authoritative adoption as well, preventing a terminal fallback report with
  no presentation owner.  Hubble shaded/wire cold/warm now
  pass on both System GL and OSMesa; retain this as the root-overview contract.
- **P1 user feature:** `view cutting` owns a disabled-by-default world-space
  plane (`enable`, `origin`, and `normal`) separately from `view zclip`; the
  controller and exact-pick clip list are covered by CTest.  The qged View
  Settings faceplate now exposes the enable toggle and world-space
  origin/normal fields through the public GED command, and its focused test
  proves widget-to-command and command-to-widget synchronization, preserves a
  command-updated plane across unrelated faceplate changes, and rejects a
  zero normal without a partial plane update.  Plane readback is intentionally
  canonical (`normal * distance`), so tangential components of a point are
  not retained as a misleading second state.
  **Implemented presentation 2026-08-24:** an enabled plane now retains a
  camera-projected, post-CAD HUD grid.  Its extent follows the synchronized
  camera, it never owns clipping state, and it is deliberately post-CAD so the
  plane cannot clip or occlude its own cue.  The focused controller test
  verifies enable/disable and a subsequent camera synchronization; System GL
  and OSMesa captures confirm the grid itself is visible on both backends.
  Normal navigation must never enable it.

  **Resolved retained-render parity 2026-08-24:** a compiled System-GL CAD
  assembly formerly ignored Coin's active clipping state because its Tier-2
  shader/indirect paths bypass fixed-function clip planes.  While one or more
  Coin clip planes are active it now uses the retained fixed-function VBO
  route, which consumes the same accumulated `SoClipPlaneElement` state; the
  normal Tier-2 path returns immediately after sectioning is disabled.  A
  scripted moss replay verified an off/Tier-2 frame, a visibly partial
  enabled/Tier-1 section, and an enabled/then-disabled return to Tier 2.
  The equivalent OSMesa capture presents the same retained half-space.

  This is an intentional correctness fallback for the optional feature, not
  the final high-performance clip implementation.  Before adding direct plane
  manipulation or claiming large-section performance, pass explicit
  application-neutral world-plane data through `CadViewState` to every shader
  tier, the CPU wire path, and exact rendering, then add a backend-tolerant
  pixel regression for transformed/instanced geometry.  The affordance remains
  a presentation/picking client of `view cutting`, never a second clipping
  authority.
- The significance frontier now orders optional transitions by projected error
  reduction, footprint/silhouette significance, user emphasis, complete
  transition cost, and the current frame budget.  `libBObol_lod_update_action`
  exercises one prominent occurrence against sixteen cheap ordinary
  occurrences, proving that the sole affordable richer prefix goes to the
  prominent one.  Still qualify this on constrained real large scenes:
  prominent wheels, blades, tails, and hulls must not remain conspicuously
  coarse merely because many cheap objects exist.
- **Implemented terminal admission refinement 2026-08-25:** complete source
  discovery now publishes an immutable profile (unique asset/occurrence
  counts, encoded source-byte total, largest asset, and reuse evidence)
  through the compact stream to the live source.  The allocator uses it only
  to extend its bounded perceptual planning window from 2,048 to 4,096
  occurrences for a complete, low-cardinality, bounded-byte source.  Unknown,
  large, or large-asset sources retain the incremental planning path.  The
  stream contract test covers profile handoff; `ged_test_drawing_lod` and
  `ged_test_obol_draw_sync` pass with the consumer enabled.

  This is a conservative gate, not a substitute for the view allocator.
  Only a visible, materially sized occurrence may use the direct
  full-mesh fast path, and it must still compete with the same protected
  visual-importance and draw/memory budgets as PoP.  Subpixel occurrences
  aggregate regardless of source size.  A 150k visible scene must therefore
  remain on compact contracts and bounded importance-ranked PoP/aggregate
  admission rather than gaining a source-size-only direct escape hatch.  The
  terminal branch now performs that full decision: the controller enables it
  only for a complete safe scene profile; the action rechecks exact visibility,
  aggregate-point eligibility, BoT source type, supported draw mode, and the
  remaining live render-cost allowance.  Work enters the ordinary bounded LoD
  service.  The worker publishes immutable prepared CAD geometry, repeated
  occurrences share that geometry and one assembly part, and display/result
  byte accounting makes terminal assets subject to ordinary memory eviction.
  PoP preparation may continue in the background because later large-view or
  memory-constrained operation still needs it; it is no longer a prerequisite
  for the first terminal presentation of an affordable scene.

  Focused coordinator, submit-action, and full-detail-provider tests cover
  rejection, direct provider selection, immutable worker preparation,
  retained reuse without resubmission, and fallback.  The real qged matrix
  passed cold and warm shaded
  Generic Twin runs on both System GL and OSMesa; each backend reached complete
  visible-target coverage without lingering boxes.  Large profiles remain
  excluded by construction and must continue to pass the 50k/150k matrix.
  `ObolLodAdmission.tla` models these modes and passed TLC on 2026-08-25
  (52 generated / 40 distinct states, depth 5).  Its first run caught the
  otherwise plausible error of admitting a cheap ordinary direct mesh before
  an affordable prominent mesh; preserve protected-floor ordering when this
  reaches the implementation.
- The bounded protocol proof is in `ObolLodArbitration.tla`: it covers coverage
  priority, affordable prominent-floor non-starvation, budget safety,
  unchanged-epoch monotonicity, and explicit terminal constrained debt.  Run
  TLC with `ObolLodArbitration.cfg` after changing the arbitration protocol.
  The current model passed on 2026-08-23 (102 generated / 72 distinct states,
  depth 12), together with the allocation oracle and LoD update/coordinator
  tests.
- `test_bobol_retained_allocation_oracle` exhaustively checks the bounded
  numeric marginal frontier and deterministic minimax tie handling against an
  independent oracle.  Keep it with the TLA+ model when changing allocator
  policy.  Images remain the authority for perceptual quality.
- `ObolProgressivePipeline.tla` replaces the former cross-policy mode model.
  It treats few-small, many-small, single-large, distinct multi-large, and
  shared-large scenes as profiles of one pipeline rather than control modes.
  TLC passed on 2026-08-25 (329,996 generated / 191,250 distinct states, depth
  37), proving bounded preparation ownership, evidence-gated planning,
  truthful terminal outcomes, active progress witnesses, eventual termination
  after finite input, and non-starvation of distinct large assets.  It does
  not prove real-time FPS or image quality; retain perf and GUI evidence.
- `ObolLiveSpatialPublication.tla` is the complementary producer/publication
  guard for cold spatial-cache construction.  It checks that a completed page
  is the only page eligible for publication, coverage remains until the final
  hierarchy, cancellation freezes the publication set, and a durable cache
  marker cannot escape a cancelled or incomplete hierarchy.  TLC passed it on
  2026-08-24 (343 generated / 212 distinct states, depth 11).  This makes the
  TLA+ suite a useful check on composition of behaviour modes, but it does not
  predict wall-clock cost, frame rate, or perceptual mesh quality; retain the
  allocator oracle, perf, and image-matrix qualification for those concerns.
- Preserve the resident-allocation boundary.  Quiet residency stops at the
  allocator-selected presentation cut; active scale interaction may prefetch
  one bounded transition past it.  Reintroducing unconstrained quiet prefetch
  recreates the 150k allocation loop; disabling it during scale interaction
  strands Lucy.
- Preserve the closed-population transaction: do not allocate a scene-wide
  population while provider, task, result, publication, residency-drain, or
  interrupted-render work can still change it.  A terminal capacity witness
  must not rearm without a view, policy, population, or capacity edge.
- **Resolved focused liveness counterexample 2026-08-25:** result batching and
  frame acknowledgement now share one typed presentation transaction, while
  the HUD is an observer of derived convergence state.  The exact 16 GiB-
  capped 150k cold shaded replay passes on System GL and OSMesa in 26.4 and
  17.8 seconds: all 150,001 leaves are available, zero boxes or failures remain,
  producer/result/submission/frame work is empty, and the terminal outcome is
  explicitly CONSTRAINED.  This closes the known crash/liveness counterexample;
  it does not replace cold/warm wire or visual-significance release rows.
- Establish explicit release thresholds for first useful image, interaction
  latency, static convergence, stable-frame time, peak/resident memory, and
  visual-quality-floor debt.  Record p50/p95/max rather than one favorable run.

## P0: bounded cold source preparation

- **Resolved source-admission undercount 2026-08-24:** the outer realization
  coordinator previously charged every root at a fixed 256 MiB before a
  worker entered the database.  That was unsafe for multiple large direct
  leaves: each could immediately deserialize far more data than its admission
  reservation.  A zero caller estimate now resolves from the immutable
  snapshot path and directory before worker admission.  Leaves reserve their
  serialized record size times the bounded import-copy factor plus fixed
  overhead; a combination reserves the coordinator's finite capacity because
  its own directory record says nothing about descendant imports.  Explicit
  estimates remain available for specialized callers and tests.  The source
  realization contract test creates a real disk leaf and proves automatic
  admission reserves the derived bytes before its worker callback starts.
  This prevents the outer multi-root undercount; it does not replace the
  separate byte-bounded page/cache-builder contract below.

- **Investigated 2026-08-24:** a fresh System-GL Lucy GUI-matrix cold run
  held its cache-producing PoP task in flight past the ordinary 45-second
  `wait_progressive_idle` deadline and accumulated a 958 MiB incomplete
  cache.  A verified empty-cache timing replay establishes the actual
  foreground progression: coverage bounds at 0.4 s, global 28.0M-face / 14.0M-
  point PoP classification at 3.8 s, spatial chunks at 11.8 s, and a visible
  340k-face whole-object payload at 25 s.  The GUI test must distinguish that
  useful retained payload from ongoing durable cache construction; do not
  classify it as an empty-view failure or hide a true cache nontermination by
  merely extending every interaction deadline.  The cache remains a P0 memory
  and completion concern: establish an explicit cache-growth/admission bound,
  then requalify cold/warm Lucy.  Earlier evidence is
  `/tmp/qged-lucy-system-bounded-20260824` and
  `/tmp/qged-lucy-cold-timing-cold.json`; incomplete cache directories were
  deliberately removed.

- **Resolved harness/endpoint liveness 2026-08-24:** the GUI driver now has
  `wait_progressive_payload_ready`, which waits for a retained CAD payload
  without treating optional cache persistence as an input barrier.  A true
  cold, pre-created `BU_DIR_CACHE` Lucy replay produced its first useful mesh
  at 8.25 s (559,494 faces, one payload, zero structural boxes) and became
  terminally idle at 57.7 s.  The corresponding warm replay produced its
  first mesh at 0.28 s and became idle at 1.76 s.  The cold terminal transition
  exposed a System-GL render-only latch: a no-pump final HUD/presentation
  request could remain coalesced in a nested event loop.  `QgCanvasState` now
  presents that narrow terminal state synchronously, while active progressive
  work remains queued and deadline-bounded.  Once all visible targets have
  current-quality witnesses, an in-flight reusable-cache build is now a
  `BACKGROUND` convergence state rather than false `REFINING`: the HUD says
  "Building reusable LoD cache  view usable", reports task count, and shows an
  indeterminate active segment.  Only a terminal, non-background `IDLE`
  snapshot with a complete fraction may say "View ready"; the observation
  boundary falls back to "Finalizing view" if those fields ever disagree.
  This proves cache isolation and payload-vs-
  idle behavior; it does **not** discharge the required explicit cache-growth/
  admission bound or complete Lucy matrix qualification.

- Cold source preparation must be admitted before it calls `rt_db_get_internal`.
  Under a 12 GiB address-space cap, System GL Lucy reached that librt path and
  `bu_calloc` invoked `bu_bomb`; the existing working-set policy permits an
  exceptional oversized task to run alone, which is not a safe guarantee.
  Replace that exception with a source-preparation admission contract that can
  return a constrained, diagnosable result without entering a non-recoverable
  librt allocation path.  It must preserve a useful existing presentation and
  permit retry after a capacity edge.
  **Resolution 2026-08-24:** detached database-source realization now refuses
  an estimate above its coordinator limit before it queues a worker.  The
  terminal `CONSTRAINED` result preserves structural coverage and is explicitly
  logged by GED's deferred-realization owner.  The source-realization contract
  proves no callback or worker starts for an over-limit request.  The shared
  mesh-provider governor now has the same rule: an over-limit task is removed
  from the queue only to publish a terminal diagnostic, without invoking its
  provider or charging an impossible active reservation.  The LoD-service
  contract covers both the process and per-service rejection paths.  This
  closes the unsafe oversize bypass; bounded source-page construction and
  durable-cache completion remain separate work below.
- **Reconfirmed 2026-08-24:** a focused empty-cache Lucy System-GL probe
  published a 559,494-face whole-object payload with zero structural boxes in
  1.16 seconds, but the one active cold task reported a 7,178,315,648-byte
  peak estimated working set while the transient governor was configured for
  1 GiB.  That was the former single-oversize-task bypass in
  `lod_task_working_set_available` and `bobol_lod_working_set_acquire`, not a
  measurement anomaly.  It is now rejected before provider execution; rerun
  the constrained Lucy row to establish the user-visible structural-coverage
  and retry behavior.  Do not use the old quick-preview result as current
  cold-source qualification.  Bounded importable/chunked source preparation
  and durable-cache completion remain required for useful constrained detail.
- **Spatial-producer integration gap 2026-08-24:** the serialized spatial
  producer is cache-tested and can keep its reservation within the governor,
  but its first 4,096-face source-order seed page is not a whole-object mesh.  A
  qged trial correctly reported zero boxes and a bounded 1 GiB reservation,
  yet the rendered image was only a tiny isolated Lucy fragment.  Do not
  enable that route as normal production policy until its early publication
  is coverage-aware: retain the overview/leaf coverage until a global
  representative mesh exists, or construct a bounded globally representative
  spatial preview.  The cold task must then reserve its bounded producer
  working set, followed by durable completion, cold/warm zoom/turnover, and
  constrained-address-space qualification.
- The first guard is implemented and verified: the formerly crashing 12 GiB
  System GL Lucy shaded/wire replays now terminate without `bu_bomb`, with one
  constrained occurrence retaining its structural fallback.  Address-space
  refusal is reported as terminal `FALLBACK`, rather than an implementation
  error, so it is diagnosable without poisoning provider-failure accounting.
  It is safety evidence only.  Cold large-mesh preparation still needs a chunked/importable
  source representation so Lucy can provide a useful mesh rather than that
  fallback under this constraint.
- A V5 serialized-source bridge now handles the first safe alternative when
  raw import admission fails.  It borrows the memory-backed database record,
  validates its envelope and indices, converts points in fixed bounded
  blocks, and uses recoverable C++ allocation failure rather than entering
  `rt_db_get_internal`/`bu_bomb`.  It preserves authored winding/culling and
  deliberately declines BoTs with active authored surface normals until their
  optional records are included in the direct-source contract.  The cold
  serialized path is regression-tested after clearing the cache on a reopened
  memory-backed V5 database.  This is a safety/availability bridge, **not**
  completion of the paged source-preparation P0: native point/face arrays and
  the existing global PoP activation structures are still materialized for one
  asset.  Qualify the 12 GiB Lucy route before claiming it improves that case,
  then replace those global arrays with the byte-budgeted page contract below.
- Do not mistake `BObolStagedSourceMesh` for that representation.  The compact
  discovery detail worker deliberately skips authored BoTs, and its staged
  lease is a short handoff window; the later LoD task consequently reimports
  Lucy through `rt_db_get_internal`.  Extend the existing compact producer
  with a bounded authored-BoT page/import path, or persist an equivalent
  coarse PoP prefix before that lease ends.  It must remain view-demanded and
  must not eagerly import every authored BoT during hierarchy discovery.
  The shared V5 serialized-BoT reader now safely exposes fixed vertex/face
  records to the recoverable bridge, but the cache builder still converts them
  into whole native arrays.  The missing seam is therefore a streaming
  cache-builder/source-page contract with an explicit byte budget, not another
  borrowed-view fallback or a longer staged-source lease.
- The page contract must be deliberately narrow and measurable:
  1. a source reader supplies validated indexed point/face records from either
     native arrays or serialized V5 bytes; its native fast path must remain
     pointer-based;
  2. a serialized reader keeps a bounded per-worker decoded vertex/face page
     working set, with no process-lifetime thread-local pages;
  3. cache identity must remain the canonical post-winding point/face stream,
     so a safe source path cannot create a visually identical but redundant
     hierarchy merely because it arrived from serialized bytes;
  4. per-corner authored normals either participate in that reader or decline
     cleanly—never silently alter shaded output; and
  5. the final source lease and all pages are released before quiet residency
     management begins.  Initial source paging need not remove the separate
     global activation/cache-output structures; account for those explicitly
     in admission and reduce them in a subsequent builder pass.
- **Implemented 2026-08-23:** the V5 source reader is now used through PoP
  classification, initial-prefix materialization, spatial chunk construction,
  and cache serialization.  A serialized producer retains one 16,384-record
  point page and face page per active phase/worker; native sources remain
  direct pointer loads.  The test clears the cache, cold-generates from the
  reopened memory-backed V5 record, proves cache-key equality with the native
  path, and checks CW culling.  Active authored surface normals and
  degenerate serialized faces still decline cleanly rather than changing
  shaded semantics.  This removes full native source-point/face arrays from
  the constrained path, but not the global activation/cache-output structures.
- Acceptance for that page contract is a cold V5 Lucy run under the historical
  12 GiB address cap with a mesh preview (not only a box), no `bu_bomb`, a
  bounded measured page working set, stable post-generation source residency,
  and byte-identical cache identity versus the ordinary native source path.
- **Next producer design — spatial PoP leaves:** retain the existing consumer
  contract of one hierarchy with globally stable quantization cuts and
  independently drawable chunks, but stop deriving those chunks from a
  source-wide activation graph.  After the bounded global vertex-bounds scan,
  partition source faces deterministically into fixed spatial leaves.  A leaf
  owns its source-face page list, local vertex remap, per-cut face/point
  counts, bounds, and immutable payload record.  It may be classified and
  persisted independently under a byte reservation; no leaf build may retain
  another leaf's decoded points, indices, or remap.  The hierarchy descriptor
  is published only after its leaf table is complete, while each leaf record
  has its own completeness marker.  View admission may then request visible
  leaves immediately and schedule their refinement independently.

  The first implementation must keep one global quantization schedule and
  use source coordinates at shared boundaries, so neighboring leaves render
  the same snapped position at a given cut.  It must not invent transition
  triangles or an independent per-leaf coordinate system in the first pass.
  That gives crack-free equal-cut rendering and lets the current renderer,
  which already consumes chunk-local records, remain intact.  Neighboring
  unequal-cut stitching and adaptive meshlet error metrics are a subsequent
  enhancement, not a reason to retain global source arrays now.

  Required invariants and tests:
  1. cache identity is the canonical post-winding source stream plus spatial
     partition/version, independent of producer-worker timing;
  2. every valid source face occurs in exactly one leaf; local indices are
     in range and each leaf's terminal record is exact for its face set;
  3. aggregate terminal face count equals the sanitized source count; shared
     vertices may be duplicated across leaves and therefore only have a
     lower-bound aggregate point-count invariant;
  4. a leaf is independently recoverable after interrupted publication; no
     partial leaf is discoverable from a hierarchy descriptor;
  5. cold scheduling can publish a bounded visible leaf mesh without waiting
     for distant leaves, while a warm hierarchy remains view-selective; and
  6. 12 GiB Lucy must prove these facts through System GL visual capture,
     with a bounded peak and no process-wide termination.
- **Producer constraint confirmed 2026-08-23:** chunk serialization no longer
  allocates source-vertex-sized remap tables; each existing page now keeps only
  a leaf-local source-id map.  The page-record writer and its metadata writer
  are now separately callable from a local source-face/activation stream; the
  existing global producer uses that same path and the mesh-cache suite passes.
  Do not enable chunk-only cache publication atop the current
  global builder: a 20 GiB Lucy GUI replay died immediately after global chunk
  construction when page-record writes overlapped the still-live global
  activation/prefix structures.  Page records must instead be emitted by the
  bounded leaf producer, after the prior leaf's scratch is released.  This is
  a peak-memory invariant, not a cache-format preference.
- **Implemented test-gated 2026-08-23:** the serialized source can now run a
  deterministic spatial producer (`BOBOL_LOD_SPATIAL_LEAVES=1`).  It streams
  fixed-width source-face records into anonymous fixed-grid spools, drains one
  bounded page at a time through the local page writer, and publishes only the
  final metadata/marker after all page records succeed.  Its unit fixture
  verifies both a durable completed route and an injected constrained route:
  the completed descriptor has a partition-specific cache key, zero global
  cluster ranges, exact terminal face coverage, independently readable pages,
  survives reopening, and publishes one bounded direct first-page preview;
  the constrained route has nonzero live coverage but remains undiscoverable
  after its generating handle is released.  `ObolSpatialProducer.tla` models
  this producer boundary: direct preview while a cache write is held, atomic
  durable discovery, process-local constrained pages, cancellation, and no
  cache re-entry during a write transaction.  TLC passed its two-leaf bounded
  model on 2026-08-23 (35 generated / 28 distinct states, depth 6), including
  the pre-identity mapped-source bootstrap lifecycle.  It is
  deliberately not the
  default yet: cache identity/versioning, partial live publication after a
  capacity edge, and real Lucy cold/warm qualification remain before it may
  replace the existing useful-prefix path.
  The producer now emits a 4,096-face direct first-page handoff before the
  complete face stream is spooled; durable leaves remain at the ordinary
  65,536-face target.  A System-GL Lucy diagnostic replay reached that handoff
  in 16.1 ms of producer time and the service accepted it (`published=1`),
  eliminating the former cache-lock self-deadlock.  It still missed the
  replay's 1.5-second image checkpoint because task/scene scheduling reached
  the handoff after that checkpoint, and it did not settle inside the
  intentionally short five-second deadline while full classification
  continued.  Treat this as evidence that first-page construction is bounded,
  not as first-useful-image qualification.  Follow-up tracing separated two
  gates: the controller's coverage-only census correctly protects a
  multi-object scene, while serialized generation was hashing and
  bounds-scanning the full source before it could materialize a first page.
  A test-gated sole-source exception now permits the spatial request while
  the census continues, and a nonterminal 4,096-face sampled bootstrap is
  published before global cache identity/bounds preparation.  The 2026-08-23
  System-GL Lucy replay confirms that the service accepts it (`published=1`,
  4,096 scene faces).  An initial whole-model-stride bootstrap took 1.95 s
  because it decoded almost every mapped vertex page; changing it to one
  contiguous source page reduced that work to 9.8 ms.  The remaining
  first-useful-image bottleneck was the roughly 1.3 s structural discovery
  pass, not page caching, mapped point/face decoding, or GUI adoption.  The
  serialized AABB scan is now a bounded parallel reduction only for a large
  individual BoT (at most two scans, four workers each); ordinary many-leaf
  discovery retains its existing outer worker pool.  The repeated real Lucy
  System-GL replay reduced coverage to 845 ms and recorded the accepted
  4,096-face bootstrap by the 0.9-second checkpoint.  Re-qualify cold/warm
  50k and 150k cases before making a broad discovery-performance claim; add a
  visible spatial-frontier representation before enabling the route by
  default.
  A 2026-08-23 20 GiB qged run selected ordinary native preparation, while
  the corresponding 12 GiB process terminated before LoD preparation began;
  neither is evidence for or against the bounded producer.  Do not attribute
  that pre-provider process failure to page classification.  The disposable
  GUI cache bundles were removed after recording this fact.
- **Implemented test-gated constrained publication 2026-08-23:** a page write
  failure can retain a strictly capped 64 MiB immutable live-page set for the
  generated handle.  No name mapping, descriptor, or completeness marker is
  published in that state, so reopening remains a cache miss; the current
  handle can still read its retained pages through the ordinary chunk callback.
  The mesh-cache test injects the write failure and proves nonzero drawable
  coverage, source-limited status, and that the incomplete hierarchy is not
  discoverable after the handle is released.  This is the correct cache
  atomicity boundary.  It still needs view-prioritized first-page selection and
  a structural representation for unbuilt pages before it provides the final
  cold Lucy user experience.
- **Investigated 2026-08-23:** the serialized reader is correct enough to
  classify the full Lucy source under a 20 GiB diagnostic cap: 28,040,274
  faces and 14,020,180 points completed bounds, classification, and attempted
  persistence in 15.9 seconds.  The same route still fails the 12 GiB
  acceptance case before a mesh is published.  This is not an occurrence-route
  rejection or a malformed V5 record.  The current constrained LMDB map retry
  is 512 MiB; a complete Lucy hierarchy exceeds it, and the monolithic cache
  transaction treats persistence failure as source-generation failure.  Do
  not paper over this by raising the virtual map under the address cap or by
  weakening the cache-completeness marker.  Design a bounded cache-publication
  protocol which keeps incomplete records unreachable, allows an already
  materialized minimum prefix to be published as a constrained live asset,
  and records exactly which higher prefixes are unavailable for retry after a
  cache-capacity edge.  Its compatibility with existing header readers must
  be proved by the mesh-cache suite before it replaces the current atomic
  transaction.
- **Implemented but not yet qualified 2026-08-23:** immutable cache records
  now commit independently and the final marker remains the completeness
  witness.  If the durable publication fails, the producer retains its
  view-bounded initial prefix as a live asset, caps that asset's advertised
  hierarchy at the retained cut, and reports it as memory-limited instead of
  publishing a partial name mapping or returning to boxes.  Focused
  mesh-cache and LoD-service tests pass.  Removing a redundant serialized
  face-validation sweep allowed the 20 GiB System GL replay to classify Lucy
  in 4.3 seconds, construct spatial chunks in 13.1 seconds, reject the
  undersized durable cache safely, and visibly retain a detailed stable mesh
  at `/tmp/qged-lucy-paged-fastpath-20g-20260823`.  The 12 GiB replay still
  terminates before classifier completion, so this is **not** the required
  acceptance evidence: reduce the producer's remaining whole-mesh working
  set or add a demand-driven spatial-leaf builder before claiming constrained
  Lucy is fixed.  Keep the source-limited asset contract, but do not use it
  to mask that earlier failure.
- Do not treat a cache-map fallback as the solution.  The cache now has normal,
  constrained, and minimum mapping retries, but its capacity affects reuse,
  not safe raw-BoT materialization.
- **Resolved cache-map admission seam 2026-08-24:** the mesh-cache map is now
  clamped to one quarter of a finite Linux `RLIMIT_AS`, including an explicit
  `BOBOL_MESH_LOD_CACHE_GB` override.  An unavailable bounded map remains an
  ordinary cache miss rather than falling through to libbu's default map size
  and consuming the process address-space headroom reserved for the database,
  source pages, and renderer.  `libBObol_mesh_lod_cache` passes both normally
  and under a 1 GiB process cap with an 8 GiB requested cache map.  This
  constrains cache address-space admission; it does not change the separate
  cold-source working-set or durable-cache-completion requirements above.
- **Corrected 2026-08-24:** the bounded serialized spatial producer remains
  explicit opt-in through `BOBOL_LOD_SPATIAL_LEAVES`; the service carries the
  choice through `BObolMeshLodPreviewRequest`, but does not select it from a
  million-face threshold.  A prior note incorrectly described that policy as
  automatic.  This restraint is intentional: a 4,096-face contiguous source
  bootstrap reached System GL in 1.64 seconds and OSMesa in 3.01 seconds with
  no structural fallback boxes, but it is a sparse isolated fragment in an
  all-model Lucy view.  It is not acceptable to replace useful whole-object
  coverage with that preview merely to improve a timestamp.  Automatic use
  remains P0 until a bounded, globally representative preview is defined and
  qualified for cold/warm, constrained-memory, and view-prioritized behavior.
- **Reinvestigated 2026-08-24:** the service had selected the serialized
  entry point for the explicit flag without copying that decision into
  `BObolMeshLodPreviewRequest`, so `BObolPopState` silently took its costly
  non-spatial path.  The request now propagates the flag and a short System GL
  Lucy replay reaches the 4,096-face callback (`published=1`).  The capture
  still consists of sparse surface specks, not a useful whole-object mesh.
  Increasing the sampled page to 65,536 faces took 22.9 seconds before the
  bootstrap could publish.  Both alternatives are rejected as a production
  answer.  The remaining design needs a distinct, coverage-aware preview
  representation which remains additive to the standing overview until it has
  a visual-coverage certificate; it must not masquerade as a one-cut PoP
  hierarchy.
- **Progress 2026-08-24:** the cache callback now carries an explicit preview
  kind.  A serialized source publishes a bounded, spatially stratified
  coverage sample; libBObol turns that sample into an immutable 24-cells-per-
  axis occupancy-surface preview for the current shaded or wire channel.
  It is a direct coverage payload, not `BObolLodProgressiveMesh`.  The former
  partial spatial-leaf callback had no consumer and was removed rather than
  retaining a redundant page-sized mesh copy.  The service supplies already-discovered source
  bounds so preview construction avoids a duplicate pre-preview bounds scan;
  durable generation still independently validates its stored bounds.  Fresh
  System-GL Lucy captures reached the coherent preview at approximately 3.36
  seconds with 26,410 bounded source representatives; shaded showed the
  coarse global silhouette and `draw -m0` showed its occupancy grid.  The
  sampler now uses at most eight workers after the source exceeds one million
  vertices, with bounded roughly-32-MiB sampling scratch; a fresh OSMesa Lucy
  run constructed coverage in 1.10 seconds and displayed it at the three-second
  checkpoint.  The cache test now uses a unique per-process cache directory,
  and two consecutive runs plus the focused cache/service pair pass.  This is
  encouraging P0 progress,
  not release qualification: qualify System GL/OSMesa cold and warm behavior,
  multi-object pressure, source-bound mismatch rejection, visual usefulness
  at several orientations, and terminal handoff before enabling the spatial
  producer automatically.
- **Current System GL probe 2026-08-24:** with a newly created cache and the
  explicit spatial-producer request, Lucy published 26,410 coverage points in
  1.08 seconds.  Its three-second shaded framebuffer had a coherent global
  occupancy silhouette, one active CAD/mesh payload, and zero visible
  structural boxes; the report correctly described the cache write as one
  background task.  The intentionally short driver did not exit after its
  seven-second checkpoint because it waits for that durable task, so it was
  stopped after a bounded 67-second observation and its 236-KiB temporary
  cache was removed.  This is a presentation probe only, not a terminal or
  cache-growth result.  In particular, it demonstrates why the durable
  spatial producer must become independently bounded and cancellable rather
  than making a useful preview wait for full cache persistence.
- **Cancellation seam hardened 2026-08-24:**
  `BObolMeshLodPreviewRequest` now carries an optional cooperative
  cancellation callback.  The serialized producer polls it at bounded
  vertex/face/page intervals, including parallel global classification;
  `BObolLodService` binds it to the owning generation.  Cancellation therefore
  stops an obsolete source build without asynchronous thread termination and
  without committing a cache-completeness marker.  The mesh-cache regression
  cancels during serialized cold work and verifies both that the producer
  returns no hierarchy and that no cache payload is discoverable.  This makes
  view replacement/teardown safer, but it does not make a single spatial page
  cheaper or discharge the explicit durable cache-growth bound.
- **Regression coverage 2026-08-24:** `libBObol_mesh_lod_cache` now feeds a
  deliberately stale discovery extent into the serialized spatial producer
  and verifies that the published coverage bounds are re-derived from source
  data.  The initial OSMesa `draw -m0` replay produced the same coherent
  occupancy-wire silhouette but only at 7.39 seconds.  The bounded read-only
  sampler now partitions only serialized vertex sampling across at most eight
  workers (about 32 MiB maximum sample scratch, plus bounded reader pages); it reduced the measured
  coverage pass from 2.90 s to 1.22 s and put that same preview at the 3 s
  checkpoint.  Fresh System GL and OSMesa shaded probes showed the solid
  silhouette there; fresh System GL and OSMesa `draw -m0` probes showed the
  matching occupancy wireframe.  An OSMesa probe with the resident limit set
  to 64 MiB had one active coverage payload and `view_ready` by 3.34 seconds.
  After removing the unused leaf callback and raising only the bounded sampler
  cap, a fresh OSMesa shaded probe measured 1.10 s for coverage sampling with
  eight workers and again showed the silhouette at three seconds.
  A simple 20-second OSMesa hold retained that payload at both checkpoints.
  A foreground exact-current cold lifecycle invocation was terminated by the
  test runner after its `first-payload-ready` checkpoint without a
  report/result row; repeat it under a managed long-lived runner before
  treating cold/warm handoff as qualified.
  This is a real improvement, but not release clearance:
  validate shaded, warm-cache, constrained-memory, multi-object, and System
  GL behavior before automatically selecting the producer.
- **Rechecked 2026-08-24:** that early handoff is not presently a
  user-visible whole-object answer.  A fresh System-GL Lucy cold lifecycle in
  `/tmp/qged-lucy-profile-contract-20260824` retained the structural overview
  at 0.33 s and 1.94 s, first published a CAD payload at 59.4 s, and reached
  a quiet 2.10M-face presentation at 61.2 s.  Its completed 1.30 GiB cache
  then produced a mesh at 0.34 s and a quiet 2.10M-face presentation at
  2.84 s in `/tmp/qged-lucy-warm-audit-20260824`.  This is correct
  structural behavior--the source-order 4k page remains hidden rather than
  falsely replacing the whole model--but it confirms the cold globally
  representative publication is still P0.  Do not use the producer callback
  timestamp as first-useful-frame evidence.
- **Live-page transport progress 2026-08-24:** the cache producer now has a
  typed `BObolMeshLodSpatialPageCallback`, distinct from the whole-source
  coverage callback.  It emits only a fully validated deterministic page at a
  common terminal cut after its cache record is complete.  The cache test
  verifies cumulative local indices, deterministic page order, full source
  face coverage, and cancellation after one delivered page with no durable
  hierarchy marker.  `BObolLodPresentationLayer` is now carried and validated
  by `BObolLodResult` and `CadPayload`; duplicate or invalid layer identities
  are rejected.  This is producer/payload infrastructure, not live page
  rendering: the compact adapter still has one active part per logical
  occurrence.  It must gain bounded auxiliary renderer-instance ownership
  before this callback is connected to `BObolLodService`; do not publish a
  local page as the sole occurrence geometry.
- Requalify System GL Generic Twin, Lucy, Hubble, and Havoc after this path is
  hardened.  Generic Twin passed cold/warm shaded/wire at 12 GiB; Hubble wire
  materialized payloads but its first-frame checkpoint was empty, while Hubble
  shaded retained a render/refinement-frame liveness witness.  Neither Hubble
  result is release evidence.

  **Progress 2026-08-23:** the Hubble OSMesa shaded cold replay initially
  exposed a test-host liveness gap after zoom: all service queues were empty,
  but the nested `wait_progressive_idle` loop observed an outstanding render
  request without asking the visible endpoint to consume it.  The event driver
  now presents such requests through the same real-canvas path it uses for
  scope barriers; it does not synthesize controller completion.  The rebuilt
  cold replay passed in 17 seconds, including the Hubble prominent-object
  floor assertions.  A separately launched certified-warm replay then passed
  in 19 seconds.  It exposed one validator inconsistency: Hubble's selected
  panel may be depth-occluded, and the comments already allow a zero erase
  delta, but the redraw threshold was still nonzero.  The semantic
  compact-entry/selected-instance checks prove the exact erase/redraw
  transition, so the Hubble redraw threshold now matches the documented
  occlusion rule.  The current System GL shaded cold/warm pair subsequently
  passed in 9/8 seconds with the same camera-contract and prominence-floor
  checks (`/tmp/qged-hubble-system-20260823`).  The matching System GL wire
  pair and OSMesa wire pair also passed; Hubble is now qualified for this
  cold/warm, shaded/wire, camera, selection, and subpath-redraw matrix on both
  renderers.  Other models and the complete System GL matrix remain open.

## P0: release qualification

- **2026-08-24 current-binary evidence:** `libBObol_lod_update_action`,
  `libBObol_mesh_lod_cache`, `ged_test_drawing_lod`, and
  `ged_test_drawing_select` pass.  The selection test now keeps unpacked APNG
  controls in `ged_select_test_controls` and uses the sibling
  `ged_select_test_runtime_cache` for mutable facetization/LoD state; do not
  put controls beneath `BU_DIR_CACHE`, because a valid cache clear otherwise
  destroys later image comparisons.  Hubble shaded System GL cold/warm passed
  all camera, zoom, selection, and subpath erase/redraw checkpoints in
  `/tmp/qged-hubble-cutting-requal-20260824`.  The equivalent OSMesa pair
  passed with a 15-second per-view qualification allowance in
  `/tmp/qged-hubble-cutting-osmesa-15s-20260824`; Hubble's many-component
  software path can exceed a diagnostic 3-second window but did converge,
  with no pending work, at every checkpoint.  Generic Twin's exact-current
  full compact-path matrix is now green in
  `/tmp/qged-generic-current-20260824`: all System GL and OSMesa cold/warm
  shaded/wire rows plus their camera-contract checks passed (11--21 seconds
  per graphical run).  This supersedes the partial
  `/tmp/qged-generic-requal-20260824` evidence, but does not qualify the
  remaining real-model, resize, or constrained-memory rows.
  **Rechecked after cold-preview cancellation hardening 2026-08-24:** a fresh
  true-cold/certified-warm shaded Generic Twin row passed the full deterministic
  qged event sequence and camera comparator on System GL in 16/14 seconds and
  OSMesa in 34/34 seconds, respectively.  The temporary evidence directories
  were `/tmp/qged-generic-system-current-20260824` and
  `/tmp/qged-generic-osmesa-current-20260824`; their essential result rows are
  retained here and the disposable images/cache logs may be removed.
  Hubble shaded cold/warm is now current-binary green on both renderers in
  `/tmp/qged-hubble-current-20260824` (System GL 9/8 seconds; OSMesa 19/14
  seconds), and its wireframe counterpart is green in
  `/tmp/qged-hubble-wire-current-20260824` (System GL 9/8; OSMesa 17/16
  seconds).  The matrix's prominent-payload, quality-floor, camera, and
  terminal-idle checks passed; the OSMesa terminal capture reports its
  explicit responsiveness-limited state rather than cycling.
  **Rechecked System GL shaded 2026-08-24:** the exact-current Hubble
  cold/warm row again passed the deterministic hierarchy, selection, and
  camera sequence in 16/13 seconds.  Its temporary capture directory was
  `/tmp/qged-hubble-system-current-20260824`; retain this result row rather
  than the disposable screenshot/cache files.
  Its LoD-auto resize row also passed on both renderers in
  `/tmp/qged-hubble-resize-current-20260824` (6/12 seconds), covering active
  and stable resize, minimize/maximize/fullscreen restoration, a resize storm,
  LoD toggle round-trip, and tree selection before/after the storm.  Broader
  hierarchy/selection qualification remains open.
  Havoc's exact-current cold/warm shaded/wire matrix passed on both renderers
  in `/tmp/qged-havoc-current-20260824` (7--11 seconds per run), including
  every camera-contract check.  This closes its basic renderer/mode row;
  resize, LoD-off, evaluated-mode, and edit/selection lifecycle coverage are
  still separate requirements.
  The qged GUI matrix now has its previously missing `nist` case, using an
  explicit `BOBOL_NIST_DB` override or an installed `NIST_MBE_PMI_*.g`
  fixture.  Its current cold/warm shaded/wire System GL/OSMesa run passed in
  `/tmp/qged-nist-current-20260824` (5--12 seconds per row), with terminal
  tessellated BREP payloads, no structural boxes, and no pending work.  This
  is basic BREP renderer/mode coverage; adaptive growth/zoom-out reclamation
  and a large real BREP remain separate P1 qualification.

- **Cold large-single-mesh preview remains open:** the serialized spatial
  producer is now selected explicitly in `BObolMeshLodPreviewRequest`, not
  by a process environment variable, and has durable/constrained cache unit
  coverage.  It is not normal production policy yet.  A 4k contiguous source
  page remains a producer-local diagnostic only: it is never eligible to
  replace a whole-object presentation.  The current typed coverage callback
  instead emits a bounded spatial sample, which the service converts to an
  immutable occupancy-surface payload.  That payload is globally bounded by
  the source extent and has passed short cold Lucy shaded/wire probes on both
  renderers.  Automatic selection still requires full lifecycle, memory,
  turnover, and visual-utility qualification; never regress to promoting a
  local source page merely to improve a first-image timestamp.
  The completed implementation is specified in
  `live_spatial_publication.md`: it publishes immutable completed pages as
  assembly layers beside coverage, never replaces coverage with a local page,
  and does not rebuild an aggregate mesh on the GUI thread.  Warm cache reads
  and cold production now converge on the same retained-page representation.
  `ObolLiveSpatialPublication.tla` now makes the required ordering
  independently checkable (343 generated / 212 distinct states on
  2026-08-24): page completion precedes publication, coverage persists until
  final complete geometry, and cancellation cannot publish a cache marker.
  This is a control-plane guard for the implemented publication contract, not
  evidence of renderer throughput or perceptual quality.
  **Progress 2026-08-24:** a focused true-cold trace found the actual
  publication blocker.  Lucy completed global PoP classification in about
  4.2 seconds, but its complete 559,494-face / 1.20M-point prefix was
  rejected against the provisional 500k first-publication allowance; the old
  path therefore waited for chunk construction and cache persistence.  The
  sole-source bootstrap now permits one byte-capped global preview, and the
  service selects the richest *contained global cut* which fits the current
  renderer's scene allowance.  It narrows the producer's borrowed prefix view
  without copying or rebuilding it.  System GL published the 559,494-face cut
  with zero structural boxes.  OSMesa, whose measured allowance was 12,352,
  published the globally ordered 5,656-face / 12,884-point cut at 8,941 cost,
  with a 4.15 ms frame and zero structural boxes.  The latter is deliberately
  coarse but preserves the whole Lucy silhouette; it is not the unsafe local
  source page.  The GUI matrix now requires a cold Lucy first payload within
  15 seconds.  Re-run the complete cold/warm lifecycle after this change; a
  short trace is evidence for the first-publication contract only, not final
  convergence, turnover, or cache-memory qualification.
  **Qualified System GL shaded 2026-08-24:** the exact-current cold/warm
  lifecycle at `/tmp/qged-lucy-system-current-20260824-rerun` passed every
  matrix checkpoint and camera-contract check (92 s cold, 36 s warm).  The
  cold first payload was the globally ordered 5,656-face cut with zero
  structural boxes; the terminal cold and warm views both had one retained
  payload, 7,410,570 active faces, zero structural boxes or pending tasks,
  a 1.0 convergence fraction, 230,334,788 resident bytes, and a 0.241-pixel
  maximum certified projected error.  This discharges the current System GL
  shaded first-publication/lifecycle regression.  It does not replace the
  required OSMesa, wire, constrained-memory, real-model, or platform rows.
  **Qualified OSMesa shaded 2026-08-24:** the matching exact-current
  cold/warm lifecycle at `/tmp/qged-lucy-osmesa-current-20260824-rerun` also
  passed every matrix and camera-contract checkpoint (93 s cold, 43 s warm).
  Its cold first payload was the same globally ordered 5,656-face cut with no
  structural boxes.  Both terminal views had one payload, 340,082 active
  faces, no pending task or structural box, a 1.0 convergence fraction,
  13,698,584 resident bytes, and a 1.265-pixel error certificate.  That
  lower terminal cut is the declared software-renderer quality/performance
  tradeoff, not an unresolved presentation state.  Wire, constrained-memory,
  and broad real-model qualification remain open.
  **Qualified System GL wire 2026-08-24:** the cold/warm row at
  `/tmp/qged-lucy-system-wire-current-20260824` passed all lifecycle and
  camera-contract checkpoints (93 s cold, 37 s warm).  It replaced the
  overview with a global 5,656-face cold preview and settled both runs at
  7,410,570 faces, zero structural boxes, no pending task, and a 0.241-pixel
  certificate.  OSMesa wire and constrained-memory coverage remain open.
  **Qualified OSMesa wire 2026-08-24:** the initial cold row found a real
  return-history gap: a complete 106 ms static frame narrowly missed the
  nominal 100 ms deadline, so it was not remembered and `autoview` rebuilt
  rather than recalling its existing quality.  History now accepts a complete
  static image within a bounded 20 percent scheduling tolerance; it remains
  only a seed and interactive frames retain their strict deadline and normal
  re-admission.  `libBObol_lod_coordinator` covers the zero, normal, and
  overflow deadline cases.  The true-cold/warm rerun at
  `/tmp/qged-lucy-osmesa-wire-current-20260824-fixed` passed all lifecycle
  and camera-contract checkpoints (115 s cold, 52 s warm): 5,656 cold-preview
  faces, 225,112 terminal faces, zero boxes/pending work, a 1.585-pixel
  certificate, and one exact history recall.  Constrained-memory and broad
  real-model qualification remain open.
  `ObolColdPresentation.tla` is the bounded policy model for this seam.  It
  proves that a local bootstrap cannot hide the overview, a mesh publication
  requires global coverage, and either a mesh or an explicit constrained
  overview eventually terminates visually; cache completion is background
  work after that point.  TLC checked both reservation outcomes on 2026-08-24
  (7 generated / 7 distinct states, depth 4).  It cannot establish a
  wall-clock deadline or judge whether a representative preview looks useful,
  so retain the real Lucy captures as the acceptance authority.

- **Resolved 2026-08-24 regression:** Generic Twin hidden-line (`draw -m4`)
  over OSMesa now satisfies the LoD resize convergence contract.  The
  completed root-cause record is in the visual-quality section above: after a
  point classifier had settled, broad eligibility still caused the controller
  to abort the same irreducible software frame.  The terminal row passes with
  zero boxes and no pending handoff.  Do not reintroduce an eligibility-only
  deadline veto; deadline work must correspond to an actionable population
  transition.

- **Resolved 2026-08-26:** multi-occurrence capacity samples no longer consume
  their bounded invalid-sample allowance while a renderer-wide handoff ceiling
  is still hiding the allocation being calibrated.  The applied occurrence
  plan is now sufficient to complete that handoff even when the allocating
  pass changed cuts; the exact ceiling-free successor frame owns the sample.
  A ceiling at or above every active cut is classified as inert rather than as
  a different population.  Focused coordinator tests cover changed-plan
  handoff, the applied-but-hidden candidate transition, and exact capacity-
  population classification.  Before the final transition fix, four identical
  warm OSMesa multi-Lucy/xpush replays split evenly between pixel demand and a
  premature ordinary-deadline endpoint.  After the fix, six of six early-
  checkpoint replays reached the identical pixel-demand allocation in
  4.59--7.04 seconds: approximately 984,100 faces, 1,548,889 render-cost units,
  0.457-pixel maximum projected error, no presented structural boxes, and no
  terminal owner or obligation.  The eight progressive payloads expand to
  4,230 presented spatial occurrences; payload count is therefore not a
  triangle or draw-count proxy.  The corresponding System GL run reached
  3,807,054 faces at 0.228 pixels in 2.45 seconds.  This fixes the prior
  timing-dependent coarse endpoint; cold eight-asset preparation and the full
  turnover/memory matrix remain open.

- **Resolved 2026-08-26:** the current format-22 warm OSMesa Lucy path no
  longer reproduces the obsolete rectangular first-candidate slab.  The full
  196-sample qged lifecycle at `/tmp/qged-current-qualification` presents a
  recognizable 9,336-face whole silhouette at the 0.2-second checkpoint,
  reaches coherent detail by 1.5 seconds, and passes zoom, rotation, lighting,
  hierarchy selection, subpath erase/redraw, camera, quality, and terminal
  contracts.  It finishes with 2.10M faces, 0.632-pixel certified error, zero
  boxes, and no owner or obligation.  Retain the visual oracle: this result
  closes only the warm shaded row; cold and wire qualification remain open.

- Run the exact-current graphical matrix from truly cold and certified warm
  caches.  The exact-current 50k OSMesa shaded/wire rows are green; the 150k
  OSMesa shaded/wire cold/warm rows are green under the bounded
  diagnostic environment.  Re-run them on production hardware scale rows
  before release.
  **Resolved regression 2026-08-26:** the current Generic Twin shaded matrix
  passes cold and warm on System GL and OSMesa at
  `/tmp/qged-generic-cross-renderer-current-20260826`.  Exact-view return
  recalls history and improves rather than degrades the maximum projected
  error on both renderers; every terminal checkpoint is ready, box-free, and
  occurrence-failure-free.  This retires the 2026-08-24 coarse-return result
  as historical evidence rather than current debt.
  Re-run representative real-model rows if their cache or policy path changes.
- Run the complete System GL matrix on a usable X host.  Do not promote older
  System GL captures to current-binary evidence.  Compare semantics and images
  with OSMesa within declared tolerances; use APNG and apitrace for a flash,
  corruption, camera jump, or GL-state mismatch.
- Qualify Generic Twin, Lucy, multi-Lucy/xpush, Hubble, Havoc, NIST BREPs,
  Stanford meshes, and independent multi-gigabyte vehicle models.  Exercise
  shaded/wire/evaluated modes as applicable, LoD on/off, resize, zoom,
  rotate, translate, selection, exact/subpath erase-redraw, and cold/warm
  cache behavior.
- Run ASan/UBSan across the shared dynamic stack with active worker teardown,
  cache corruption, endpoint replacement, plugin reload, edit cancellation,
  and rapid view close.  Run TSan/LSan on a compatible native worker.  A
  sandbox runtime limitation is not a successful race/leak check.

## P0: interaction and editing qualification

- Drive actual Qt controls and mouse paths, not only equivalent GED commands.
  For every interaction, widget state, GED semantic state, retained scene,
  and command readback must agree after GUI- and command-originated changes.
- **Validated 2026-08-23:** the dedicated fractional-DPR qged selection replay
  passed on both OSMesa and System GL using the real view-selection palette.
  It covers point select/remove, held drag-rectangle placement in physical
  canvas coordinates, selection-list readback, and selected erase/redraw.
  The runner now accepts either the installed or regression `moss.g` fixture,
  so a space-conscious partial build does not silently lose this coverage.
  This is one matrix row, not clearance for hierarchy-scale selection or the
  remaining controls below.
- **Rechecked 2026-08-24:** the focused interaction gate passed 22/22 current
  binaries: librt descriptor conformance; compact/edit-manipulator ownership;
  GED selection, faceplate, and transaction semantics; Qt primitive preview,
  selection, measurement, and faceplate synchronization; MGED restore; and
  actual qged event, polygon, sketch, primitive-edit, framebuffer, and
  settings replays.  This is implementation-substrate evidence, not the
  required real-renderer, fractional-DPR, plugin-lifecycle, or hierarchy-scale
  completion matrix.
- **Rechecked 2026-08-24 after the current LoD policy work:** all 35
  `edit_runtime` handlers and the 15 focused GED/Qt/qged/MGED interaction
  checks passed.  The latter includes selection semantics, primitive preview,
  faceplate/selection/measurement synchronization, qged event replay,
  polygon and polygon-sketch replay in single and quad layouts, primitive-edit
  replay in single and quad layouts, settings, and MGED edit restore.  This
  confirms the current implementation gate; it does not replace actual mouse,
  fractional-DPR, lifecycle, or hierarchy-scale qualification.
- **Rechecked 2026-08-24:** all seven registered qged event replays pass on
  the current build: event-script dispatch, polygon creation/manipulation in
  single and quad layouts, polygon/sketch synchronization in both layouts,
  and primitive-edit interaction in both layouts.  This is deterministic
  in-process GUI coverage, not evidence for physical-pointer, fractional-DPR,
  plugin-lifecycle, or hierarchy-scale behavior.
- **Resolved command-readback race 2026-08-24:** qged copied a command result
  only after view/database observer delivery.  A faceplate refresh can itself
  execute `view cutting`, so a successful `view polygon area` command could
  report that unrelated multi-line status instead of its numeric result.  The
  model now snapshots `ged_result_str` immediately after `ged_exec`, before
  observer delivery.  The qged polygon single/quad replays, along with the
  focused 17-test GED/Qt/qged/MGED interaction subset, now pass.  Preserve
  this command-result boundary for all observer-driven GUI controls.
- **Rechecked MGED Obol baseline 2026-08-24:** shaded-mode, GQA prefix,
  host-configuration, smoke, ERT, progressive-LoD, and edit-restore rows pass
  (7/7 available tests).  The Tk/Obol attach and framebuffer rows are
  capability-skipped in this headless environment; they remain explicit
  graphical-host qualification work and are not counted as passing coverage.
- Complete the qged selection, polygon, measurement, and faceplate matrix:
  point/rectangle selection and modifiers; tree styling and draw/erase under
  selection; polygon creation/manipulation/booleans; all measurement modes;
  lighting/navigation/axes/grid/framebuffer controls; resize and fractional
  DPR in single and quad layouts.
- Complete the primitive/sketch work described in `qged_editing.md`: the
  librt descriptor contract now classifies all registered operations, but
  rejection/readback behavior, the reusable manipulator vocabulary,
  specialized sketch curves, and common qged/MGED/gsh sessions including
  sed/oed still require qualification.

## P1: scale, memory, and maintainability

- Measure discovery on cold local storage, warm page cache, and slower storage.
  Keep parallel read-only discovery bounded by I/O throughput and memory, and
  prove ordinary librt/GED hierarchy and tree publication parity.
- Qualify repeated-instance visibility turnover separately from distinct-asset
  stress.  Verify spatial-page reuse, zoom-out compaction, and cache reload
  without boxes or full asset rebuilds.
  **Current evidence 2026-08-26:** the complete shared multi-Lucy System GL
  replay passes with eight occurrence payloads, one asset producer, zero final
  boxes, and stable close-view/rotation turnover.  It exposed and now guards
  the distinction between requested spatial pages and prepared renderer pages;
  an off-view or unprepared occurrence can no longer drive an endless
  structural-repair loop.  The xpushed eight-asset cold row is making bounded
  progress with four workers under the 1 GiB aggregate working-set limit, but
  exceeds the current 150-second first-cache deadline.  Two interrupted runs
  grew the resumable cache from 942 MiB to 1.8 GiB and promoted an additional
  asset without OOM.  Warm initial-view System GL and OSMesa capacity rows now
  pass with exact occurrence allocations as recorded in the 2026-08-26 visual
  quality entry above.  Complete and profile the cold cache, then run the
  System/OSMesa visibility-turnover and memory-recovery scripts before closing
  this row; do not reinterpret an incomplete cold run as a pass.
- Qualify adaptive BREP tessellation/LoD growth and zoom-out reclamation over
  the NIST corpus and at least one large real BREP model.
- Keep extraction ownership-led.  The first extraction pass exists, but
  `database_source.cpp` (about 19k lines) and `draw_obol.cpp` (about 13k) are
  still concentrated.  Extract only at a demonstrated lifecycle or profiling
  seam; preserve compact storage, sparse deltas, bounded owner-thread windows,
  and immutable payload ownership.
- Replace the remaining monolithic static `ged_test_obol_draw_sync` coverage
  with owning-library tests plus a production shared-stack GED integration
  test.  Static linkage belongs only in explicit static-link checks.

## Exit condition

The stack is production-ready only when the exact final binaries pass the
release matrix on required renderers and platforms, no visual correctness or
responsiveness blocker remains, real large-model quality meets the agreed
significance floor, memory recovers after view turnover, and qged/MGED/gsh
editing and interaction remain coherent under lifecycle stress.
