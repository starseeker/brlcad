# QGED object and primitive editing

Last reviewed: 2026-08-24

This document defines the remaining production work for qged object, primitive,
and sketch editing.  It is not a compatibility plan for the early primitive-
specific plugins.  Qt controls, commands, MGED interaction, and retained Obol
manipulators must all operate on one renderer-neutral GED/librt edit session.

## Authority and ownership

```text
librt descriptor + rt_edit implementation
                    |
                    v
       path-scoped GED edit session
          |          |                 |
          v          v                 v
       edit CLI   Qt controls    retained manipulators
          ^          ^                 ^
          +----------+-----------------+
                 typed events
```

- Librt owns primitive semantics, validation, edit modes, keypoints, feature
  topology, operation descriptors, current-value readback, and interaction
  math.
- Libged owns path resolution, occurrence matrices, edit-scope authority,
  revisions, transactions, checkpoint/revert, commit/cancel, database events,
  selection reconciliation, and observer delivery.
- LibBObol/Obol own retained previews and manipulators, constant-pixel
  presentation, picking, capture, and view/model motion conversion.  A node
  carries stable operation and feature identity, never mutable librt structs.
- Libqtcad/qged own generated controls, optional custom adapters, input routing,
  and readback presentation.  They never keep an authoritative primitive copy.
- MGED, gsh, Qt widgets, commands, and scene manipulators are peers over the
  same session; none reconciles a private edit state after the fact.

## Session invariants

- One canonical occurrence path has at most one writable intermediate state per
  GED context.  An alternate occurrence of the same leaf conflicts explicitly.
- A successful update increments a monotonic revision after the state is valid.
  A rejected operation changes neither primitive, dirty state, readback,
  revision, nor observer stream.
- Begin, update, checkpoint, revert, commit, cancel, conflict, and invalidation
  are typed events.  Reentrant delivery is queued or coalesced, not recursive.
- Observers reread the session after an event; event payloads do not constitute
  a second state representation.
- Commit produces one database mutation transaction.  Cancel leaves the
  database unchanged.  A clean commit does not publish a false modification.
- Database close, replacement, rename, remove, nested-path erase, selection
  changes, view replacement, and plugin unload cannot leave stale previews,
  manipulators, captures, promoted draw scopes, or borrowed endpoints.
- Editing an already drawn occurrence reuses its retained source and LoD
  residency.  Promotion must not rebuild or copy unrelated large meshes.
- All Coin mutation and picking happens on the view owner thread.  Worker-made
  immutable geometry is published atomically.
- Command, widget, label, preview, manipulator, selection, and every registered
  view report the same target, operation, values, revision, and terminal state.

## Implemented foundation

The following are substrate, not remaining design alternatives:

- installed libged edit sessions, typed submission and observers, revisions,
  preview snapshots, transactions, and command integration;
- canonical operation names, descriptor validation, typed parameters,
  capability metadata, dynamic topology bounds, and generated-value readback;
- generated Qt editor/session controller and shared libged primitive preview;
- migrated ELL, extrude, revolve, BoT, ARB, and sketch adapters without direct
  qged database writes or private editable primitive copies;
- ELL retained handles, scalable indexed ARB and BoT selection/manipulation,
  and a shared multi-view retained geometry path;
- libqtcad `QgSketchEdit`, which submits through GED and handles general
  non-unit sketch planes without borrowing `rt_edit`;
- semantic edit invalidation, rename following, alternate-occurrence conflict,
  multi-context isolation, and synchronous plugin lifetime barriers; and
- command/widget/manipulator replay coverage for the migrated adapters under
  System GL and OSMesa in single- and four-view layouts; and
- explicit coordinate metadata for every declared ARB, BoT, and sketch
  interaction.  Descriptor validation rejects inferred coordinate space so a
  retained manipulator cannot silently acquire ambiguous drag semantics.

These protections must survive subsequent cleanup.  In particular, do not
restore per-plugin primitive ownership, command-string submission from Qt,
whole-mesh rollback on pointer motion, per-element Coin nodes, or duplicated
previews.

The current focused CTest gate exercises qged polygon and primitive-edit UI
replay, qtcad preview, selection, picking, measurement, faceplate, and
multi-view synchronization.  Primitive replay projects live retained handles;
ELL, ARB, and sketch use the endpoint input path, while the BoT vertex path
uses ordinary Qt mouse press/move/release events.  That is substrate evidence,
not full interactive qualification: a current System GL host, fractional-DPR
resize, real mouse coverage for every retained manipulator family, and the
broader operation/rejection matrix below still remain.

The focused gate was rechecked on 2026-08-24: 22/22 tests passed, including
the single- and quad-view polygon/sketch and primitive-edit replays, direct
qged synthesized-input paths, selection, faceplate, measurement, framebuffer
controls, and MGED edit restoration.  It proves that the covered
command/widget/retained-scene paths still agree; it does not replace the
renderer, DPI, operation-classification, or real-user interaction
qualification listed below.

The immediately relevant focused subset was rerun after the composed LoD
policy model was added: all 35 `edit_runtime` suites and the 10 targeted qged,
Obol, and MGED edit/selection/polygon/measurement regressions passed in 19.05
seconds total.  This is regression evidence only; it does not retire the
real-mouse, fractional-DPR, or lifecycle qualification work below.

## Remaining primitive contract work

### Operation classification and behavior audit

The descriptor-classification baseline is complete.  On 2026-08-23,
`rt_edit_descriptor_contract` iterated all 34 registered descriptors and
validated 109 generated, 24 action, 89 custom, and 5 explicitly unsupported
commands.  It requires canonical names, generated-operation readback, and an
explicit coordinate space whenever an interaction descriptor exists.  This is
the owning-librt gate; qged must not duplicate it with primitive-name tables.

Every installed librt edit operation resolves to exactly one class:

- **generated**: declarative parameters and complete current-value readback;
- **action**: typed input which intentionally has no persistent current value;
- **custom**: a specialized adapter is required but state remains observable;
- **unsupported**: intentionally unavailable with a diagnostic.

No operation with an interaction descriptor may retain `INFER` coordinate
metadata.  Continue to populate selection domain, coordinate space, parameter
semantic, dynamic legal bounds, and manipulator hints in librt as new
operations are exposed.  Do not create primitive-name or topology policy
tables in qged.

The remaining work is behavioral: verify that every advertised handler's
descriptor still matches its runtime legality, readback, and retained control
path, particularly for transitional primitive handlers.

The existing per-primitive runtime suites are now registered as the
`edit_runtime` CTest label instead of being build-only executables.  The gate
contains ARBN, datum, annotation, revolve, torus, ellipsoid, TGC family,
EPA/EHY/ETO/HYP/RPC/RHC/particle/REC/superell/cline, extrude, metaball, HRT,
ARS, pipe, DSP, VOL, EBM, halfspace, sphere, ARB8, BoT, sketch, B-spline,
BRep, NMG, combination, and generic string/current-value coverage.  REC is
also included in the generated command-ID header, so the handler, test, and
clients share its named operations.  On 2026-08-24 the gate passed 35/35 in
0.35 seconds.  This makes the established valid-operation and selected
error-path tests a CI contract; it does not claim the exhaustive invalid-input
matrix below is complete.

### Rejection and readback audit

For every advertised handler, prove:

1. valid set followed by get returns the canonical value;
2. invalid scalar, vector, index, list, topology, or file input is rejected;
3. rejection leaves geometry, mode, selection, dirty state, and revision
   unchanged; and
4. descriptor bounds and handler legality agree for every primitive subtype.

Audit `ft_set_edit_mode` side effects, especially EBM, DSP, extrude, CLINE,
PIPE, NMG, and other transitional handlers.  Mode selection and operation
execution must not be doubled.  Fix the owning librt handler; a generic
whole-primitive clone/rollback is too costly for interactive large BoTs.

## Retained manipulator completion

Build a small instanced vocabulary rather than primitive-specific scene trees:

- point and indexed point sets;
- constrained axis arrows and radius handles;
- plane handles;
- rotation rings;
- indexed edge and face sets; and
- topology/custom surfaces where an aggregate picking structure is required.

Each handle has stable identity, hover/active/disabled state, accessible color,
constant-pixel sizing, depth policy, and an explicit interaction constraint.
Picking returns semantic feature identity.  The controller alone converts a
drag to typed session input.  Capture must survive resize and view replacement,
and must terminate cleanly on command-side cancel, target invalidation, tool
handoff, and plugin unload.

The implementation must scale without one QObject or Coin node per mesh
element.  Large BoT pointer motion patches only affected shared points and
does not copy or reconstruct topology.

## Sketch and profile completion

- Add direct viewport construction for lines, arcs, Beziers, NURBs, open and
  closed profiles.
- Add specialized arc/NURB affordances, control-point/weight/knot interaction,
  segment splitting, joining, deletion, and topology validation.
- Keep stable source-segment identity when curves are tessellated for display.
- Explicitly arbitrate input ownership between the 2-D sketch surface, 3-D
  object/solid manipulators, camera controls, selection, and polygon tools.
- Exercise skewed, non-unit sketch planes, unit changes, snapping, holes and
  islands, persistence, cancellation, and database invalidation.
- Replace extrude/revolve's transitional profile helper with the shared sketch
  adapter.  A linked profile edit must update every dependent retained preview
  without creating another editable copy.

Once custom lifecycles no longer require separate binaries, collapse the
migrated primitive plugins into one registered edit tool/factory.

## QGED interaction qualification

Tests must drive actual Qt controls and synthesized mouse/keyboard input; an
equivalent GED command is useful coverage but is not proof of GUI viability.
Run System GL and OSMesa, single and quad layouts, DPR 1 and fractional DPR,
and resize during an active dirty session.

For each representative generated and custom primitive:

- enter from a sole scene/tree selection and from an explicit nested path;
- change state from CLI, widget, and retained manipulator in alternating order;
- assert command `status/get`, control values, labels, preview, handle state,
  every view, and database state at every revision;
- cover create, select, move, resize, rotate, scale, and topology operations;
- cover checkpoint/revert, commit, cancel, undo, tool handoff, target switch,
  view switch, clear/zap, resize, and database close/reopen;
- cover repeated instances, alternate-occurrence conflicts, erase-under-edit,
  source replacement, rename/remove, selection/highlight changes, and
  selection falling from one occurrence to zero or many;
- prove invalid input is immutable and produces a useful UI diagnostic; and
- unload/reload the real plugin with an active session, hidden panel, visible
  modeless dialog, retained surface, and pointer capture.

Use Generic Twin, Hubble, havoc, Stanford meshes, NIST BREPs, and mixed
50k/150k hierarchies to prove editing does not restart unrelated LoD, copy
retained assets, or stall the GUI.  Use moss/rook for evaluated CSG and legacy
editing modes.

## MGED and command-line qualification

- Put `sed`, `oed`, the `edit` command, and gsh on the common session semantics.
- **Rechecked 2026-08-24:** the current MGED regression subset passes for
  `sed`, `oed`, primitive edit, rejection, knob, rotate, scale, and translate
  edit paths (8/8).  Together with the Obol edit-restore row this proves the
  presently covered command surface remains intact; it does not replace the
  required nested-path, lifecycle, mouse, or gsh qualification below.
- Prove exact nested-path and subtree promotion/restoration over compact CAD
  sources without expanding the root or rebuilding LoD.
- Exercise mouse/knob motion, labels, axes and ticks, all applicable primitive
  operations, checkpoint/revert, commit/cancel, alternate occurrence conflict,
  erase under edit, and rejection images.
- Verify the exact pre-edit retained scene and selection are restored on every
  exit path.
- Retire duplicate legacy edit state only after image and interaction parity,
  not merely successful command return.

## Completion criteria

Primitive editing is production-ready only when:

- every advertised operation is explicitly classified and conformant;
- generated and custom adapters pass the same bidirectional replay contract;
- sketch curve creation and specialized topology editing are complete;
- the manipulator vocabulary covers the required semantic constraints without
  per-element scene overhead;
- qged and MGED pass real mouse interaction and image-level qualification on
  both renderers; and
- commands, controls, retained presentation, selection, database transactions,
  and all views remain one coherent state under invalidation and lifecycle
  stress.
