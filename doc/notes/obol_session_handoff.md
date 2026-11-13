# Obol drawing: session handoff

Updated 2026-09-20. Start with the [simplification guide](obol_simplification_guide.md)
for direction, [active debt](libbobol_active_debt.md) for tasks, and
[production readiness](obol_production_readiness.md) for acceptance. This is a
replaceable checkout snapshot. The [integration record](obol_main_integration.md)
owns the rebase and repair history; the [September 19 archive](obol_20260919_handoff_history.md)
preserves earlier work.

## Current checkout

The user authorized source/documentation rework toward a smaller supported
mutation surface and prioritized rebasing on current `origin/main`, including
the new `annotate` command. The rebase is complete: branch `brlobol7` uses
integrated commit `d2e8cb382c`, based on fetched main `0d745fca35`. Current
source/documentation repairs are uncommitted. No remote changes were made.

The complete pre-integration state is preserved at `123010cebd` on
`brlobol7-before-main-20260919`, including the old history and working changes.
The integration record maps the upstream behavior to its new drawing owners.
Continue repairing demonstrated defects within those owners and the guide's
bounded production workflow.

BRL-CAD uses `.build-main`. The native Obol/OSMesa build is `.build-main-native`,
installed normally into `.build-main-native/install`. Obol is the sibling
checkout at `a544ae3876bd` with a 68-path local change set, including the two new
`CadDisplayPlane` files; its OSMesa submodule is clean at `bfa90de454a7`.
The current manifest captures tracked and untracked source, build settings,
installed headers, executable/library hashes and runtime dependency resolution.
A branch name or an unpatched Obol checkout does not identify the tested stack.

Use an explicit build-local install prefix: the native build's default prefix
points to the external bundle. Never replace libraries while an application or
test is using them. Recheck loader resolution after configuration. Reproduction
commands are in the [integration record](obol_main_integration.md#reproducing-the-native-dependency-baseline).
Keep historical `.build-main-deps` and `.build` unchanged.

## Current behavior and boundaries

- Erasing an annotation retains unrelated pending source jobs. Clear/teardown
  still cancel all jobs; per-source admission rejects obsolete results.
- Direct and grouped primitives share effective material resolution. Legacy
  `rgb` and canonical `color` edits survive annotation regeneration. The color
  resolver's private absent-region sentinel does not change primitive identity.
- A late root overview cannot replace an existing leaf. An edit to a root
  overview uses complete-source publication, preserving leaf identity and
  selection; overviews cannot seed reusable primitive geometry.
- View assembly preparation derives its records from staged presentation and
  leaves canonical source records unchanged. Native indexed replacement and
  bounded sparse batches avoid a second complete occurrence vector.
- Screen annotations retain raw display offsets, model anchor and pixel density
  in immutable native geometry. One projection helper supplies rendering,
  bounds, picking, rectangle selection, collection and viewport export.
  Independent cameras can zoom, rotate and resize without changing source
  records or rebuilding the frame plan.
- GED and Qt viewport exports use `SoBRLExportAction::applyViewport`. Ordinary
  node/path export retains anchor-plus-XY layout; a completed viewport traversal
  does not leave its camera attached to a reused action.
- Live GED snapping and Qt measurement use matching viewport-scoped actions.
  Ordinary traversal remains deterministic and viewless; path-local measurement
  retains stored display offsets. The model-coordinate-only legacy vlist
  publisher returns `-1` for screen annotations, which remain on native compact
  `CadDisplayPlane` geometry.
- Filled areas retain librt triangulation, including holes, in existing native
  triangle storage with an explicit unlit role. Fills draw before strokes in
  every mode; automatic wire-mode picks include fill interiors. Source summaries
  and cache reuse retain plotted identity, and export includes fills plus lines.
- Authored widths are compact immutable native wire runs. They multiply the
  occurrence's base pixel width, survive zoom and geometry replacement, and
  reach library export as floating-point widths. Ordinary geometry avoids a
  width array; geometry identity and memory accounting include authored runs.
- The same ordered wire runs retain annotation segment color and exact 16-bit
  line pattern. Compact fill runs retain authored RGBA without changing
  triangle topology. Geometry colors render through fixed, GLSL, batched and
  direct-software paths where applicable; selection, highlighting and explicit
  draw color replace them at the occurrence. Library export retains effective
  color, fractional width, pattern mask/factor and fill color.
- `RT_ANNOT_ROLE_MASK` is retained as an immutable semantic fill flag. Fixed
  and GLSL rendering reproduce the active solid or gradient viewport background;
  selection cannot recolor it, and library/Qt export retain the mask explicitly.
- PostScript and PNG consume effective line width, exact pattern/factor,
  color/alpha and annotation fills from viewport export. Masks use page white
  or the PNG command background. Plot uses a declared named-pattern mapping and
  explicitly cannot retain widths or filled areas. Qt object queries carry
  exact style vectors parallel to line and triangle geometry.
- The finite live-mutation inventory now names the owner, lifetime, revision,
  failure/callback contract and caller class for every supported family.
  Attached source invalidation, compact selection/frontier/style and fallback
  wireframe realization publish through scene controllers. Dead GED branches
  that could attach raw source children were removed. Detached fixture builders
  remain a separate supported lifetime class.
- The convergence snapshot owns one immutable progress-display classification:
  visibility, terminal readiness and a stable publication class. Libged owns
  the faceplate's text, color and geometry; Qt owns cadence and event-loop
  scheduling. Qt no longer mirrors libged policy or caches three partial values.

The earlier raw viewport-field watcher remains removed. Detached fields and
explicit rebuild configure detached viewports; owning host/stream APIs publish
live changes. Retain the preceding publication and lifetime regressions.

## Evidence and open acceptance

API-01/S1 evidence is retained under
`.build-main/obol-qualification/20260920-api-inventory`. It records the caller
search, public/installed API checks and the focused source, selection and GED
synchronization regressions. These ownership changes do not alter annotation
pixels, so the seven current attributed controls remain appropriate and were
not replaced.

FLOW-01/S2 evidence is retained under
`.build-main/obol-qualification/20260920-flow-01`. Named headless and graphical
tests exercise cold draw, camera input while work is active, close, reopen and
redraw through production owners. The headless path observes a real delayed mesh
task. The Qt path uses a bounded service-delayed task so widget ownership is
tested independently of whether the small fixture's planner requests another
cut. Both require visible initial and changed terminal output, compact terminal
geometry, terminal convergence, empty transient queues and Linux worker
retirement. Each focused lifecycle test passed ten consecutive runs. The
combined focused/broad rows, full BRL-CAD build (1,442 Ninja steps), installed-package consumer,
both `annotate` rows, repository and license checks pass. `RESULTS.md` in that
directory records commands and limitations.

Production-library evidence is retained under
`.build-main/obol-qualification/20260919-annotation-fills`:

| Check | Current result |
|---|---|
| Native Obol CTest suites | 9 pass, including installed consumer and both software rendering backends |
| Native fill/width pixels | Fill holes, all native modes, mixed point/fill/stroke geometry, screen projection, picking and role replacement pass on both software GL backends; preceding width assertions also pass |
| Affected BRL-CAD integration rows | 39 pass, no skips |
| Focused drawing checks | 4 pass, including the screen annotation's authored-width export |
| Dedicated library/Qt export tests | 2 pass; source export checks fill-only/mixed annotations and preserves fractional widths |
| Source/state publication sweeps | Both pass; allocation and callback-order assertions retained |
| Retained mutation/view protocol models | Preceding projection candidate passes 128,904/480 distinct states; this change leaves that modeled protocol unchanged |
| Current annotation image comparisons | All 7 pass against attributed controls at the unchanged 0.99 threshold |

The fill regression originally exported 24 lines and zero area for a square
with a square hole. It now exports the expected area of 12, no compatibility
outlines or triangulation strokes, and retains separately authored lines.
Native pixels and picks verify coverage through camera and geometry changes.
The renderer reuses existing GPU triangle buffers and the unlit point pass;
direct software wire falls back to GL for parts needing this pass.

The follow-up at `20260919-annotation-framing` adds command assertions and
attributes image differences without changing production libraries. With the
old camera set through GED commands, frame seven scores 0.999907 against its
unchanged control, above 0.99. Its ordinary score remains 0.892470: the camera
change is intentional. All 2,551 truck leaves agree with librt's current material
lookup; the old lookup made them uniformly brown. New assertions protect
representative materials, authored radial spans and centered complete framing.
The [attribution record](obol_main_integration.md#annotation-image-attribution)
and diagnostic manifest distinguish these checks from full qualification.
With user authorization, frame seven's control first changed to the intended
camera contract. A second diagnostic under `20260919-annotation-framing-1-6`
then fits frames one through six to their old controls using only a camera scale
and translation. Their SSIM rises to 0.961952--0.975646; model-space annotation
color bounds match exactly or within one pixel. The remaining pixels are
independently attributed to current librt materials, fixed-pixel display-plane
projection and authored styling. All seven controls now use the stable current
render, the threshold remains 0.99, and all seven comparisons pass. The prior
APNG, diagnostic source, fitted-camera images and measurements are retained.

The post-framing focused checks add native color/pattern pixels across fixed,
GLSL, instanced, flat and direct-software wire paths, plus fill-color pixels in
all draw modes and both display coordinate modes. The libBObol prototype
retains exact line masks, fill colors, semantic background masks and selection
behavior in compact export. Native mask pixels pass with solid and gradient
backgrounds in fixed and GLSL paths, both coordinate modes and the direct
software backend. The output-consumer follow-up uses an isolated annotation to
verify PostScript width/pattern/color/fill output, plot's named dash mapping and
PNG style/fill effects; the Qt test verifies exact parallel style vectors. These
checks pass, but they are not yet a sealed replacement for the candidate-wide
evidence above. A focused interaction follow-up verifies viewport-scaled snap
and measurement, camera reset on action reuse, path-local behavior and explicit
legacy vlist rejection. The same follow-up forces immediate GL on
compatibility-profile OSMesa and System-GL contexts and passes the full authored
style pixel test at tier zero, separately from direct software wire.
The library records do not by themselves establish cross-client annotation
visual parity.

At the earlier color checkpoint, a continued seven-frame annotation run left frames one
through five and seven byte-identical to the framing candidate. Frame six now
shows previously black-on-black segments in their authored colors and improves
from 0.886535 to 0.887488 SSIM. Its control still uses different framing, so it
was not replaced. The command-only regression passes.

The subsequent mask repair leaves all seven rendered frames byte-identical to
that post-color run. It therefore does not justify another control update. One
initial standalone Xvfb invocation exited during startup; the same run completed
under GDB and on an immediate normal rerun with identical frame hashes, so no
reproducible mask-path crash is currently established.

The output-consumer repair is recorded under
`.build-main/obol-qualification/20260919-output-consumers`. The command-only and
Qt regressions pass. Its continued seven-frame run reports the same frame-one
through-six SSIM values (`0.861125`, `0.847549`, `0.846693`, `0.870093`,
`0.872103`, `0.887488`), frame seven passes exactly, and all generated frame
hashes match the mask checkpoint. That output-only repair did not itself
justify a control change. The subsequent framing diagnostic supplies the
missing independent attribution and control decision.

S0--S3 are demonstrated locally; S4--S6 remain open. Offscreen software OSMesa
is not native GPU evidence. The historical load-sensitive scene-light witness
is now closed by kind-aware semantic equality and a deterministic unused-NaN
reentry regression. Large-model/BREP, remaining sanitizer lifetime, Lucy
quality, contention, editing, client and platform obligations remain in the
active backlog and production matrix. Old temporary GUI reports and warm
caches cited in historical notes are not current qualification.

## Resume here

Begin S4 with LIFE-01 and CONTROL-01, retaining the demonstrated FLOW-01
workflow as the integration guard. OWN-01/S3 is complete for the current
inventory. Its first seam made the LoD service the owner of lock-consistent
service-wide and per-generation work snapshots. Its
complete-idle predicate covers all transient service work, and its generation
result-work predicate treats a shared-producer lease as live work even when the
consumer owns no task. GED/GUI settling and controller/renderer readiness use
those contracts. The former controller-local result-pending atomic is removed;
the service queue owns availability and the controller retains only the
first-ready timestamp used for publication batching. Focused and broader
controller/model rows, compact-publication stress, installed consumption and
ten repeated runs of both production lifecycle paths pass; evidence is under
`.build-main/obol-qualification/20260920-own-01-service-work`.

The resident-capacity candidate is also complete. Stable renderer bytes are one
atomic service accounting value updated by realization, both compaction paths,
eviction and stop. `residentCapacityStatus()` samples it with the mutex-owned
growth reservation and limit. Realization publishes stable bytes before
releasing its reservation, so an observer may conservatively see both but
cannot see neither. Renderer headroom/pressure, convergence and qged consumers
use the new contract. Arithmetic, real publication, compaction, eviction/reload
and stop checks pass; the full-hierarchy publication test observes a live
reservation, changes the limit and rejects a reservation-to-stable gap. The
focused 13 rows, broader six rows, compact-publication stress, installed
consumer and both lifecycle paths repeated ten times pass under
`.build-main/obol-qualification/20260920-own-01-resident-capacity`.

The progress-display candidate is complete. A direct classification table,
ten focused faceplate/progressive/API rows, thirteen production-focused rows,
six controller/model rows, the installed consumer and both lifecycle paths
repeated ten times pass under
`.build-main/obol-qualification/20260920-own-01-progress-display`. This repair
changes classification ownership and scheduling inputs without changing
faceplate or annotation pixels, so no control image was replaced.

The closing audit migrates GED command and GUI qualification waits to the
existing host-work snapshot while retaining service worker/cache quiescence as
a separate service boundary. Remaining synchronized view, renderer and
source-admission values have one writer and explicit transaction boundaries;
diagnostics do not drive policy. The qged build, ten focused rows, and ten
repetitions of GED cross-run plus both production lifecycle paths pass under
`.build-main/obol-qualification/20260920-own-01-s3-close`.

For LIFE-01, start with the retained source/service and borrowed host/widget
lifetime cases, then run worker-active close, endpoint replacement, plugin
reload and cancellation on the shared dynamic stack under ASan/UBSan before
the supported native TSan/LSan rows. For CONTROL-01, regenerate the quiet
planning, desktop requested-frame and load-sensitive exact-frame witnesses with
complete transition journals before changing policy. Preserve tight autoview
bounds, current librt materials, display-plane projection, viewport interaction
semantics and the shared-source boundary. A concrete new duplicate writer or
split policy reopens OWN-01; do not reopen it for ordinary boundary cooperation.

Preserve preceding evidence directories under `.build-main/obol-qualification`:
`20260919-annotation-interaction`, `20260919-annotation-framing-1-6`,
`20260919-output-consumers`, `20260919-annotation-masks`,
`20260919-annotation-fills`, `20260919-annotation-strokes`, `20260919-annotation-projection`,
`20260919-view-source-isolation`,
`20260919-annotation-materials`,
`20260919-native-baseline` and `20260919-main-integration`. The framing directory's
`RESULTS.md`, `NEXT.md` and final manifest identify the attribution checkpoint.
The last sealed historical repair is CXX-SOURCE-089 endpoint destruction under
`.build/obol-qualification/20260913-endpoint-destruction`; later raw-field work
was not sealed. Serialize heavy qualification and use fresh persistent artifacts
as described in the [GUI runner instructions](../../src/qged/tests/README.md).
