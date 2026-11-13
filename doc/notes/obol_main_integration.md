# Main integration, September 2026

Status: rebased baseline with focused build and behavioral qualification.
The registered annotation images pass; cross-client visual parity and the
broader production gates remain open.  The
results below distinguish command behavior, retained rendering, GUI checks,
and evidence that still needs investigation.
The latest [image attribution](#annotation-image-attribution) explains all
seven intentional camera/material/presentation differences and records the
authorized control updates. The focused
[output-consumer mapping](#annotation-output-consumer-mappings) records the
subsequent PostScript, PNG, plot and Qt behavior without changing those images.

## Preserved starting point

The complete pre-integration history and working changes are preserved by
`brlobol7-before-main-20260919` at `123010cebd`.  The integration branch applies
one combined drawing change to `origin/main` at `0d745fca35`, starting from the
common ancestor `6c94ae71d2`.  This incorporates 553 intervening main commits.
Integration was built and tested in a separate worktree before updating the
original branch.  The safety branch retains the complete old history and
working changes.  Nothing has been pushed.

## Behavior carried into the drawing owners

| Main change | Integration decision | Qualification needed |
| --- | --- | --- |
| New `annotate` command and automatic dimensions | Preserve creation, placement, styles, updates, and show/hide. Discover displayed inputs through drawing intents; publish database mutations through GED events. Keep existing visibility and draw appearance when updating geometry. | Command regressions, retained colors, hidden and nested updates, image comparison. |
| Annotation tessellation | Adapt the new stroke-to-mesh implementation to the branch's vlist ownership. | Annotation primitive semantics and rendering. |
| Primitive colors in shaded drawing | Resolve primitive attributes in the shared material sweep for every drawing mode. Retain combination inheritance, region-table fallback, annotation defaults, and explicit draw overrides. | Upstream color fixture against retained drawing results. |
| Attribute-preserving BREP edits | Keep the new writer and notify successful mutations through GED. Remove the superseded writer wrapper. | BREP edit and attribute regressions. |
| Transactional BoT split and optional grouping | Preserve rollback, attributes, names, and grouping; publish additions only after successful completion. | Split and group regressions, failure cleanup. |
| Decimation face provenance | Keep the updated API and attributes; retain successful output notifications. | Decimation regressions. |
| Bounded facetize workers | Keep the new subprocess and staged-transfer implementation. Retire the old worker and preserve final output notifications and batch cleanup. | Worker completion, failure, cancellation, and displayed-result refresh. |
| BREP archive serialization | Keep main's shared archive reader and its lock; retain the independent drawing tessellation lock. | Archive I/O and concurrent drawing tests. |
| Submodel transforms | Keep composed transforms and relative file resolution. Remove the duplicate branch transform function. | Submodel transform tests. |
| Memory availability and resident accounting | Keep main's current probes and tests, including Linux reclaimable-cache estimates. Remove duplicate probes introduced by the automatic merge. | Memory tests and drawing capacity behavior. |
| Voxel clipping and newer geometry algorithms | Preserve upstream bounds fixes and updated APIs. | Geometry tests and compilation. |
| MGED command interruption and GUI callbacks | Preserve worker-local search interpreters, cooperative interruption, GUI dispatch, and event servicing until workers finish. Retain drawing host teardown. | Interactive cancellation, close during work, script execution. |
| MGED edit operation and reinitialization fixes | Keep explicit operation selection and the shared persistent edit initializer, using semantic occurrence paths and the new passive view adapter. | ARB and object/solid edit regressions. |
| MGED unshared view ownership | Preserve independent view-ring copies using the display owner. | Share/unshare and saved-view tests. |
| MGED framebuffer sessions | Use the GED-owned imgstream listener and create an independent session token when opening it. Keep child-tool routing and the loopback Tcl transport. | Concurrent sessions and framebuffer round trips. |
| TclCAD interpreter isolation | Preserve interpreter-local object registries and thread-local active objects. Adapt lifecycle tests to GED-owned views and endpoints. | Independent, replacement, and GUI interpreters. |
| qged startup fixes | Keep a separate valid Qt argument vector, explicit database path, plugin suffix configuration, and startup/test commands. | CLI startup, native host and scripted GUI tests. |
| New application launcher | Use an Obol display session, imgstream menu pixels, and an endpoint input layer. Preserve manifest discovery, launching, text fallback, and pointer/key/close behavior. | Build, discovery, graphical input and resize smoke tests. |
| Framebuffer standard-stream and descriptor-zero fixes | Carry the tests into libimgstream; avoid reading or truncating output streams and keep owned file handles separate from standard streams. | Both upstream file-handle regressions. |
| Legacy display-list indexes and native backends | Retain their removal. Semantic scene indexes, retained geometry, and native host providers own the corresponding behavior. | Drawing selection, material, host ownership, and platform gates. |

The shared plot adapter now distinguishes model-space annotations, whose
anchor and plane basis are already applied by librt, from the legacy
screen-space fallback.  This removes the duplicate anchor translation that
separated dimension arrows from their geometry.  Both wire representation
paths use the same adapter.

The annotation hide/show sequence also exposed an asynchronous drawing defect:
erasing one source cancelled every pending source job, leaving unrelated labels
at their provisional bounds. Source work now survives unrelated transactions;
the existing stream admission checks cancel erased or superseded sources.
Whole-scene clear and teardown still cancel all jobs. View-local LoD invalidation
remains with scene synchronization instead of being repeated by the source
provider. The annotation regression now checks delivered strokes and hidden
source absence after hide/show and a subsequent view update. It fails on the
preceding implementation and passes with this repair.

The next annotation check exposed missing material semantics in direct
primitive realization. It now uses the same full-path material sweep as
combination leaves, including legacy `rgb` aliases and canonical `color`
precedence. Direct wire realization initializes its semantic summary before
copying it into a LoD request. Annotation writers use the canonical attribute;
regeneration accepts older aliases and retains edits made with `attr set`.
The color fixture covers direct and grouped BoT/BREP drawing, and annotation
checks cover delivered colors and legacy/canonical update round trips.

Archer plugin staging now preserves implementation subdirectories.  Its
previous flattened copies caused duplicate class definitions during GUI
initialization.  Qt canvases now deliver normalized resize events to endpoint
application layers; launcher pointer buttons use the same zero-based mapping
as Qt and Tk.

The launcher now depends on the configured Qt or TclCAD window provider for its
graphical mode.  Its former independent X11/Win32 framebuffer window and event
queue are replaced by the branch's common host and input contracts.

## Qualification results

- A fresh Debug build with Qt, Tcl/Tk, OpenGL and testing enabled builds
  libBObol, GED plugins, MGED, qged and the launcher.
- Annotation command assertions pass: help, primitive/default/inherited colors,
  automatic input discovery from drawing intents, dimension and leader updates,
  hidden/nested/erased visibility preservation, and failed-update preservation.
  `ged_test_annotate_commands` runs these independently of image comparison;
  the original seven image comparisons remain enabled and unchanged.
- Retained drawing colors pass all six mode/override combinations for BoT and
  BREP primitives.  File-backed framebuffer standard-stream and descriptor-zero
  tests pass.  Image-display and window-host tests pass, including ordinary
  framebuffer visibility, pan/zoom, resize and reopen round trips.
- Annotation primitive semantics, submodel transforms, BREP archive I/O, BoT
  split/group, memory accounting, and facetize worker/validation/recovery tests
  pass.  Facetize cleanup assertions now distinguish temporary workspaces from
  the separate persistent drawing cache.
- Nine targeted MGED regressions pass: accept, ARB, loadview, oed, primitive edit,
  reject, saveview, search-exec and sed.  Host configuration, drawing, embedded
  raytrace and edit/restore smoke tests pass.  QGED startup/event replay passes.
  The fresh original-checkout run passes 34 selected behavioral rows after
  correcting the Qt controller test's platform selection. Its software and
  window-system GL assertions now run separately; both additional controller
  and host System-GL rows pass under Xvfb, including fractional-DPR coverage.
- TclCAD independent and replacement interpreters pass both headless and under
  Tk/Xvfb, including independent view state and GUI endpoint backgrounds.
  Teardown avoids evaluating commands in an interpreter already being deleted.
- Launcher discovery, console quit, graphical startup, resize, pointer motion,
  keyboard quit and a physical click on Quit after resize were exercised under
  Xvfb.  Initial and resized images were inspected.
- State-publication and source-publication sweeps pass.  The earlier
  load-sensitive scene-light assertion had a deterministic cause: point-light
  equality inspected its unused direction payload, so an indeterminate NaN
  could force child replacement during an enablement-only update.  Kind-aware
  equality and an explicit unused-NaN reentry regression close that witness.
- The shared-stack sanitizer run found an indexed progressive-rendering
  overread in the companion Obol tree.  Its GPU upload expanded quantized
  positions through the maximum active index but left authored normals empty.
  The same path now uploads the corresponding authored normal for every
  expanded vertex; the exact BRL-CAD LoD update/render regression passes under
  ASan/UBSan.
- Repository and license checks pass in the fresh original-checkout build.
  Loader checks resolve the recorded native dependency pair without an
  environment override. Reconfiguration retains that selection.
- The first post-workflow OWN-01 consolidation makes libBObol the sole owner of
  complete service quiescence and coherent generation work observation. GED,
  qged and broad settling callers consume one lock-consistent `workStatus()`
  snapshot instead of rebuilding idle policy from independent counters.
  Controller and renderer readiness consumes `generationWorkStatus()`, whose
  result-work predicate includes shared-producer leases. The service queue also
  replaces the controller's duplicate result-pending atomic; only the
  first-ready batching timestamp remains controller-local. Thirteen focused
  rows, six controller/model rows, compact-publication stress, the installed
  consumer and both lifecycle rows repeated ten times pass; evidence is under
  `20260920-own-01-service-work`.
- The next OWN-01 consolidation makes libBObol the sole publisher of stable
  resident capacity. Stable bytes are one atomic accounting value updated by
  realization, compaction, eviction and stop. `residentCapacityStatus()`
  samples it with the service-owned growth reservation and limit; exact growth
  is published before reservation release, making any overlap conservative.
  Renderer headroom/pressure policy, convergence and qged diagnostics consume
  that contract. A full-hierarchy publication test observes a live reservation,
  changes the limit concurrently and rejects an intervening capacity gap.
  Focused, model, compact-publication, installed-consumer and repeated
  lifecycle evidence is under `20260920-own-01-resident-capacity`.
- The third OWN-01 consolidation removes Qt's independent reconstruction of
  libged progress visibility, terminal readiness and phase coalescing.
  `BObolLodConvergenceStatus::progressDisplayStatus()` now returns one immutable
  display value from the complete convergence snapshot. Libged consumes that
  value while retaining text, color and geometry formatting; Qt compares it
  while retaining cadence and scheduling. The local Qt enum, mappings and
  three-field cache plus duplicate libged predicates are deleted. Focused,
  workflow, progressive, model, installed-consumer and repeated lifecycle
  evidence is under `20260920-own-01-progress-display`. This state-only repair
  preserves annotation pixels and does not replace a control image.
- The closing S3 audit migrates GED command and GUI qualification waits to the
  immutable host-work snapshot. Service worker/cache quiescence stays a
  separate service snapshot. The remaining synchronized view, renderer and
  source-admission values have one writer and explicit transaction boundaries;
  diagnostics do not drive policy. The affected build, focused transitions
  and repeated cross-run/lifecycle checks pass under
  `20260920-own-01-s3-close`. S3 is demonstrated for this inventory.

## Initial integration limits and next acceptance

At the first integration checkpoint, all seven annotation image comparisons
differed from upstream controls. The
duplicate-translation and unrelated-job cancellation defects are repaired.
The missing model-space labels in frames four and five render again, but
framing, screen-space labels and text/fill/style rendering still differ.
The material repair restores the cyan model-space leader after hide/show and
view update. Its logs, source manifest and inspected frames are under
`.build-main/obol-qualification/20260919-annotation-materials`.
No image controls or comparison thresholds were replaced at that checkpoint.
Qualify screen-space placement, text/fill/style rendering, bounds and autoview
before claiming full annotation visual parity (QUALITY-01 and EDIT-01).

The expanded color test also exposed a cold direct-BoT draw settling with only
its overview. The priority stream delivered an updated whole-target extent
after its leaf box; both shared the root path and geometry tier. Source merge
now preserves the leaf record and its refinement request against that late
overview. A deterministic delivery-order regression fails before the repair
and passes afterward; the direct-draw test passes 20 consecutive repetitions.
The failing trace and repaired evidence are retained with the material results.

The GUI checks use Xvfb and software rendering; they do not qualify native GPU,
Windows, macOS, fractional DPR, or every client workflow.  Physical cancellation
and close during active asynchronous commands, view unsharing, and the wider
large-model matrix remain required.  Passing the sourced search-exec regression
is not evidence for all interactive interruption cases.

The original integration used a recorded existing native Obol/OSMesa pair in
`.build-main-deps`. Its qualification logs, images, conflict ledger and
dependency manifest remain under
`.build-main/obol-qualification/20260919-main-integration`. The source-built
baseline below supersedes that dependency selection. Keep the older `.build`
unchanged; its binaries do not qualify the rebased sources.

## Reproducing the native dependency baseline

The initial native source state is Obol `a544ae3876bd`, with its 44-file local patch,
and its clean OSMesa submodule at `bfa90de454a7`. The exact patch and hashes,
CMake caches, build/install logs and loader reports are retained under
`.build-main/obol-qualification/20260919-native-baseline`. A branch name or an
unpatched Obol checkout does not identify the tested source.

The subsequent view/source isolation repair extends that checkout to a 45-file
patch. Its complete patch and evidence are under
`.build-main/obol-qualification/20260919-view-source-isolation`. Use that recorded
source for that preceding candidate. Obol's indexed replacement reader lets a
per-view assembly derive records without overwriting shared source storage or
allocating another full occurrence vector. The callback-compilation regression
fails before the repair and passes afterward; all nine native suites pass.

The annotation projection candidate extends Obol to a 54-file local
patch, including the new `CadDisplayPlane.h` and `CadDisplayPlane.cpp` files.
The complete patch, including untracked files, is recorded under
`.build-main/obol-qualification/20260919-annotation-projection`. Its shared
projection is used by native rendering, bounds, picking, scene collection and
BRL-CAD rectangle selection and viewport export. GED and Qt export callers use
`SoBRLExportAction::applyViewport`; ordinary action traversal preserves the
viewless anchor-plus-XY layout. An overview-to-leaf edit uses complete source
publication, preserving authoritative identity and interaction state.

All nine native suites, four focused drawing checks, 39 affected integration
rows and both publication sweeps pass on this candidate. The green screen-facing
leader appears in frame four; frames four/five changed, while the other five are
byte-identical to the preceding candidate. All seven original image comparisons still fail. Their
SSIM values are 0.862527, 0.849728, 0.849064, 0.872230, 0.874112, 0.886741 and
0.892470 against the unchanged 0.99 threshold. This establishes the tested
projection contract, not annotation visual parity or production readiness.

Dedicated library and Qt export tests also pass. The expanded library test
checks effective primitive color against librt without changing non-region
identity, and supplies an explicit camera when picking screen annotations.
The required pinewood/havoc fixture dependencies are attached after their
sample-geometry targets exist, so building the test prepares its inputs.

The failed-provider publication oracle now compares each source with its own
pre-edit snapshot. A cached mesh may legitimately use its bare asset name where
an independently constructed reference uses a full occurrence path. Diagnostics
showed that difference existed before the edit; the edit preserved its original
state. Every snapshot field remains checked. The failed sweep and diagnostic
runs remain alongside the final results.

The subsequent annotation-strokes candidate retains authored widths as compact
immutable native wire runs. Existing occurrence style supplies the base pixel
width; authored runs multiply that value, and zoom does not scale it. Admission
validates each run, source geometry identity and retained-memory accounting
include it, and ordinary wires avoid a width array. Library line-export records
now use floating-point widths. No per-segment occurrence, scheduler or extra
presentation owner was introduced.

The first independent source regression failed with an exported width of one
instead of the authored 4.5. It now covers mixed 4.5/default/2.25 widths, default
resets, wire/shaded realization and base-width round trips through the supported
display patch API. An intermediate failure from writing a raw field after
realization was corrected in the test; live field watching was not restored.
The command regression checks the screen leader's authored width after creation
and view updates.

Nine native suites pass. Pixel checks cover fixed VBO, GLSL VBO, instanced,
flat-batch and direct software drawing, including width-only replacement under
existing part IDs in a reused context. A later diagnostic capability override
forces immediate GL on compatibility-profile OSMesa and System-GL contexts;
the same style pixel test passes at tier zero separately from direct software
wire. Four focused drawing checks, all 39 integration rows and both dedicated
export tests pass. Native GPU qualification remains open.
PostScript and Qt object-level export still need explicit segment-style handling.

The width checkpoint evidence is under
`.build-main/obol-qualification/20260919-annotation-strokes`. All seven original
image comparisons fail; frames one through six changed, and frame seven is
byte-identical to the projection candidate. SSIM is now 0.861125, 0.847549,
0.846693, 0.870093, 0.872103, 0.886535 and 0.892470 against the unchanged 0.99
threshold. Controls are byte-identical to the preceding candidate. Frames three
and four were inspected alongside the control; framing, model geometry, fills
and label layout still differ. The width repair does not establish full visual
parity. Source/state publication sweeps pass in 174.54/194.82 seconds, retaining
allocation and callback-order assertions.

From the BRL-CAD checkout, with that source in the sibling `obol` directory:

```sh
cmake -S ../obol -B .build-main-native -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/.build-main-native/install" \
  -DOBOL_BUILD_TESTS=ON -DOBOL_BUILD_EXAMPLES=OFF \
  -DOBOL_FETCH_SUBMODULES=OFF -DOBOL_BUILD_BUNDLED_OSMESA=ON \
  -DOBOL_USE_SWRAST=ON -DOBOL_USE_SYSTEM_GL=ON
cmake --build .build-main-native --parallel 6
cmake --install .build-main-native
xvfb-run -a ctest --test-dir .build-main-native --output-on-failure --parallel 1
```

The ordinary install includes both runtime libraries and the Obol package.
No copied library, patched binary, loader override or external-bundle overwrite
is needed. All nine native CTest suites pass, including the installed consumer,
unit/lifecycle tests, both rendering backends and concurrent rendering.
Reconfiguration and rebuilding leave the library bytes unchanged.

Set `BEXT_BUILD_DIR` to the prepared external-dependency build directory before
configuring BRL-CAD. It must contain its normal `install` and `noinstall` trees;
this recipe rebuilds the drawing dependencies, not that entire bundle.

```sh
cmake -S . -B .build-main -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DBRLCAD_ENABLE_QT=ON -DBRLCAD_ENABLE_TCL=ON -DBRLCAD_ENABLE_TK=ON \
  -DBUILD_TESTING=ON -DBRLCAD_EXT_DIR="$BEXT_BUILD_DIR" \
  -DObol_DIR="$PWD/.build-main-native/install/lib/cmake/Obol-2.0.0"
cmake --build .build-main --parallel 6 --target \
  qged mged brlcad-launcher ged_plugins ged_test_annotate
cmake -S . -B .build-main
```

Inspect actual loaded libraries as well as `Obol_DIR`. On Linux, `ldd` on
qged, MGED and the launcher must resolve both `libObol` and `libosmesa` from
`.build-main-native/install/lib`. The existing external bundle may also stage
an OSMesa library into `.build-main/lib`; its presence does not prove which
library the process loads. Recheck after configuration or dependency changes,
and never reinstall dependencies while a drawing process is using them.

These Linux/Xvfb results establish source/build correspondence for this local
candidate (BUILD-01/S0). On this stack, all 39 affected BRL-CAD CTests and the
drawing lifecycle regression pass after the annotation cancellation repair.
The regression's retained-stroke check fails on the preceding implementation
and passes after the repair; at that checkpoint all seven upstream image
comparisons still failed.
They do not qualify native GPU behavior, other platforms, an
installed BRL-CAD distribution, or the remaining visual and lifecycle gates.

The broader production finish line remains
[the simplification guide](obol_simplification_guide.md).  Passing this
integration's checks establishes a usable new baseline; it does not close the
remaining resource, visual-quality, scale, or platform gates by itself.

## Annotation fill representation

The annotation-fills candidate retains librt's filled triangle coverage rather
than converting each triangle to wire segments. A square-with-hole fixture
initially exported 24 lines and zero area. It now exports area 12, with neither
compatibility outlines nor triangulation strokes; a mixed fixture also retains
its separately authored line. Both wire and shaded source realizations pass.
Annotation summaries and reuse caches retain plotted identity even when their
geometry includes triangle storage.

Native `shadedIsFill` marks the existing triangle channel as an unlit filled
area. Fills share the authored-point pass and existing triangle GPU buffers;
they draw before strokes in every mode, without culling or subpixel replacement.
Automatic wire picking includes fill interiors while continuing to exclude
ordinary hidden triangle surfaces. Changing the fill role rebuilds the cached
plan. Direct software wire falls back to GL when the unlit pass is needed.
Uniform point clouds retain their bulk GPU draw; fill/immediate work keeps
bounded interruption checks.

Native pixel/picking tests cover holes, all four native draw modes, mixed
point/fill/stroke parts, model/display coordinates, zoom and same-ID role
replacement in reused contexts. The fixed and GLSL software paths pass, as do
the retained width/projection checks. Nine native suites, four focused drawing
checks, all 39 affected integration rows and both dedicated export tests pass.
Final source/state publication sweeps pass in 177.56/189.52 seconds.
The exact source/dependency/binary manifest and final evidence are under
`.build-main/obol-qualification/20260919-annotation-fills`.

All seven original annotation images and controls are byte-identical to the
width candidate, with the same failing SSIM scores and unchanged threshold.
Frame four and its control were inspected. This independent fill repair does
not close the original framing, model geometry or label-layout differences.
At this checkpoint authored fill colors/masks, remaining segment styles,
downstream output-format behavior, snapping/measurement and legacy vlist
semantics remained open. Subsequent unsealed work described by the current
handoff carries segment width/color/pattern and fill color through native
rendering and exact library export. Subsequent focused work carries semantic
fill masks through native solid/gradient rendering and explicit library/Qt
export flags. Subsequent output work declares and verifies the PostScript,
PNG, plot and Qt mappings described below. A later focused interaction repair
uses the active viewport for screen-annotation snap and measurement, preserves
viewless and path-local traversal contracts, and rejects screen annotations at
the model-coordinate-only legacy vlist publisher. Cross-client visual parity
and fill-boundary interaction semantics remain open. No production
gate is closed by either representation checkpoint alone.

## Annotation image attribution

The follow-up under
`.build-main/obol-qualification/20260919-annotation-framing` changes regression
assertions and documentation; production libraries remain the fill candidate.
It distinguishes intended behavior from defects still requiring repair.

Frame seven's original comparison scores 0.892470. Main's old `_bound_objs`
expanded each object's center by its largest full-axis extent, then fitted the
union using the largest half-axis. Applying that rule to the delivered sphere
and dimension bounds predicts the control's circle center and radius. An
isolated test executable sets that camera through ordinary GED commands: the
current renderer then scores 0.999907 against the unchanged 0.99 threshold.
Both images were inspected. This attributes that frame's mismatch to camera
policy; production autoview retains centered tight bounds and libbv's
rotation-stable sphere fit. No legacy camera mode was added.

A material audit checks all 2,551 truck leaves against librt's full-path color
resolution. They agree: 2,471 brown, 28 gray, two white, 44 dark gray and six
black. Main's former display lookup resolves every leaf to the color table's
brown, explaining why some current wires appear faint or invisible on black.
Explicit region colors take precedence over the table fallback, as in librt's
region construction. The retained root overview is outside this leaf audit.
These color differences do not establish complete truck geometry parity.

The annotation command regression now checks representative delivered model
materials, authored radius/diameter/angular spans, visible coverage and a
centered geometry extent. The existing camera and material tests retain the
lower-level contracts. An isolated old-camera variant fails the new centering
assertion, demonstrating that matching the old pixels would regress the
intended contract.

With the user's authorization, frame seven of `annotate.apng` first changed to
the inspected current-camera image. A follow-up diagnostic under
`20260919-annotation-framing-1-6` fits one isotropic camera scale and translation
per remaining frame from exact authored-color pixels. Against the old controls,
frames one through six rise to 0.967124, 0.961952, 0.962733, 0.966573,
0.970344 and 0.975646. Frames one through five reproduce model-space annotation
bounds exactly and their exact-color pixels within one pixel; frame six aligns
the magenta geometry and positions the cyan/white geometry at the old bounds.
Its reduced cyan coverage is the independently verified dashed style.

The residual pixels are the already tested full-path material colors,
fixed-pixel display-plane behavior and authored styles. The green screen leader
correctly does not scale as model geometry under the substituted old camera.
All six controls therefore changed to the stable current render. The prior APNG
is retained, extraction verifies all seven replacement hashes, the threshold
stays 0.99 and all comparisons remain enabled. The registered image,
command-only and Qt export tests pass together. Remaining cross-client visual
and interaction semantics stay in the existing annotation inventory and
production gates.

## Annotation output-consumer mappings

The focused output repair consumes the exact `SoBRLExportAction` styles without
adding another display representation. PostScript converts the cyclic 16-bit
mask and factor to dash arrays, multiplies effective width by its `-l` unit,
composites alpha over page white and emits fill triangles. PNG samples the same
mask and factor in raster order, rounds effective widths to pixels, composites
alpha over existing pixels and rasterizes fill triangles with a half-open
shared-edge rule. Background masks use PostScript white or PNG's `-c` color.
Both formats emit ordinary lines, annotation fills and annotation strokes in
that order.

Plot files cannot encode arbitrary widths or area fills. They composite alpha
over white and map the five authored patterns to their named legacy line modes;
arbitrary visible masks use `dotdashed`. Qt object queries expose exact line and
triangle styles in vectors parallel to their geometry, and global triangle
records include their effective style.

The command regression isolates a styled line and fill from the truck scene.
It verifies PostScript width, dash, color and fill operators, the plot mode,
PNG changes caused independently by line style and filled geometry, and valid
background-only cleanup. The focused Qt export test verifies style-vector order
and values. Both pass. Evidence is under
`.build-main/obol-qualification/20260919-output-consumers`.

A continued seven-frame run at that checkpoint retains the preceding scores for frames one
through six (`0.861125`, `0.847549`, `0.846693`, `0.870093`, `0.872103`,
`0.887488`) and an exact frame-seven pass. All generated hashes match the mask
checkpoint. The output repair therefore does not justify another control-image
change by itself; the later framing diagnostic above supplies that decision.

## Supported source boundary after main integration

The September 20 API audit completed the finite live-mutation inventory after
the main integration. It moved compact selection, visibility and style updates
to scene-controller owners; made source invalidation publish its revision,
stale state, owned metadata and peer frame revisions as one transaction; and
routed the exact-stamped wireframe fallback through the same realization
effects. Two unreachable GED inspection branches that could attach raw children
were deleted. Production libged, libqtcad, qged, mged and gtools callers no
longer directly write attached source fields, attach source children or invoke
the audited source mutation methods. Evidence is under
`.build-main/obol-qualification/20260920-api-inventory`.

## Complete production workflow after main integration

FLOW-01/S2 is demonstrated under
`.build-main/obol-qualification/20260920-flow-01`. The named headless workflow
draws from a fresh per-process cache, observes a real delayed mesh provider task,
changes the camera, closes a worker-active production endpoint, verifies Linux
worker retirement, reopens and redraws, and requires a changed terminal image,
compact source state and empty service queues. The named Qt workflow drives the
same sequence through a real software `QgView` and the ordinary GED observer.
It uses a bounded nonpublishing service-delayed task to isolate destruction of
the widget-owned service from small-fixture planner timing, then requires an
active terminal mesh payload and changed visible `QImage` output.

Both focused workflows pass ten consecutive runs. The full GED synchronization
row retains real denial coverage, and separate named rows retain accepted and
replaced-source deferred-result coverage. All six focused/broad rows pass after
the complete build. This is Linux offscreen OSMesa evidence; native GPU,
sanitizer, large-model and other platform obligations remain in S4--S6.

The complete BRL-CAD build also exposed two integration defects outside the
focused workflow. The librt BoT edit test used pointer syntax for a stack view
object, and the LoD coordinator test omitted its transitive public Obol target;
both are repaired. A public Obol geometry callback parameter named `emit`
collided with Qt's empty `emit` macro. The parameter is now `visit`, and the
Obol regression includes the public header with an empty `emit` definition.
The full BRL-CAD build (1,442 Ninja steps), targeted native callback test, installed BRL-CAD
package consumer, both annotation rows, repository check and license check pass.
The package consumer explicitly receives the configured installed Obol package,
matching the exported BRL-CAD package's external dependency contract.
