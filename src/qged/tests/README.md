# qged graphical validation

`qged_gui_matrix.sh` drives the real qged Qt canvas through deterministic event
streams.  It is intended for release-candidate and investigative validation,
not as a short CTest.  Each run writes screenshots, command logs, controller
timing/LoD samples, cache inventories, and a summary into a new artifact
directory.

`qged_lod_quality_matrix.sh` isolates terminal image quality.  It renders the
same scene, camera, canvas, lighting, and draw mode with automatic LoD and with
LoD disabled, verifies those comparison inputs, and records libicv SSIM and
perceptual hashes.  It also measures foreground-silhouette disagreement so a
high whole-frame score cannot hide a missing thin or visually prominent
feature.  Both exact and one-physical-pixel-tolerant silhouette disagreement
are retained: the former is regression-sensitive, while the latter separates
raster-edge/antialias variation from a genuinely missing feature.  The Big Boy
side view additionally records named front-wheel, running-gear, and boiler/cab
crops in `feature_metrics.csv`.  Scores use a documented crop which excludes
the progress bar and
caption while preserving the uncropped checkpoints as user-visible evidence.
Both LoD and LoD-off crops must contain foreground pixels; a blank image is a
failed oracle, never a perfect or acceptable comparison.
`matched_assessment.csv` classifies each view as safe-direct, unconstrained
pixel-target, or constrained.  Safe-direct views enforce their photometric
SSIM floor.  Larger views enforce exact-work, geometric-error, silhouette,
proxy, prominent-floor, and constraint-evidence contracts; shaded SSIM is an
advisory review signal because simplified normals may change illumination even
when the screen-space geometry and silhouette are correct.  A low advisory
score cannot hide a topology failure, and it cannot by itself reject a
resource-constrained but geometrically valid view.
The error contract uses `presented_max_normalized_error`, measured from the
exact cuts in the captured framebuffer.  `allocation_max_normalized_error` is
retained separately as planning evidence: after a pose-only change the
allocator may describe a coarser affordable candidate while the renderer
deliberately preserves a richer resident cut.  A planning candidate cannot
make that demonstrably richer exact presentation fail its image contract.
The runner records one SHA-256 digest for each distinct database and captures
device, inode, size, mtime, and ctime around every rendering process.  It
rehashes only if that metadata changes, and a LoD/control pair is comparable
only when the content digests match.  This avoids two extra multi-gigabyte
reads per policy while retaining the mutation guard.  Database mtime alone is
not an identity: opening a writable `.g` may update it without changing the
file's bytes, while a concurrent conversion can change the actual oracle.
Big Boy is opt-in and uses the partial `big_boy.bot` mesh hierarchy.  Its
LoD-off image is a control for allocation quality within that hierarchy, not a
completeness check against the original BREP model or failed facetizations.
Scenes whose LoD-off image would exceed the host's safe memory or frame-time
envelope use `--managed-only`.  That mode runs the identical LoD camera
sequence and writes `managed_metrics.csv` with allocation, error, constraint,
proxy, residency, and renderer-memory evidence.  Deep diagnostics also write
the bounded per-occurrence sample and worst visual-importance witnesses to
`managed_payload_metrics.csv` and `managed_outlier_metrics.csv`.  It
deliberately produces no SSIM or silhouette score: those claims require a
matched control.  Large-scene qualification combines this bounded whole-scene
evidence with full-detail controls of representative semantic subsets and
named-feature crops.
`managed_assessment.csv` checks the corresponding presentation, pixel-error,
and typed-constraint telemetry but always marks the result for visual review:
without a safe matched control it is not image-quality proof.
Every policy run records `/usr/bin/time -v` resource evidence.  On Linux,
`--memory-limit-mib` places qged in a transient user scope with `MemoryMax`
and no swap allowance; use it for intentionally hostile fixtures so a failed
memory contract cannot take down the desktop.  `--initial-view-only` is the
bounded diagnostic form for measuring one settled `ae 90 0` allocation before
testing camera-transition retention separately.

The primary matrix dimensions are:

- System OpenGL and qged's `-s` OSMesa/software path.
- A newly created, empty `BU_DIR_CACHE` with no format or data children,
  followed by a warm process using the exact same cache.  (`BU_DIR_CACHE`
  itself must exist or libbu intentionally disables caching.)
- Shaded and wireframe PoP presentations.
- A deterministic 1100×800 top-level resize through the real Qt layout and GL
  resize path.  Replay mode uses a process-private QSettings scope, so neither
  saved window geometry nor plugin enablement can alter the run, and the run
  cannot overwrite the operator's settings on exit.  Camera-critical resize
  events request a bounded `stable_ms` interval; generation-tagged retries
  prevent delayed native configure acknowledgements from restoring an older
  scripted size.  The initial barrier requires a full second of native geometry
  stability, and the resize matrix records the geometry again at the semantic
  draw command so a late configure cannot silently change autoview.
- Initial draw return, 50 ms, 200 ms, 1.5 seconds, and stable checkpoints.
- `ae 90 0`, zoom in/out/return, sustained held-button rotation checkpoints,
  and post-release stable checkpoints.
- Atomic command transactions for compound camera setup.  In particular,
  orientation plus `autoview` is committed before Qt updates are re-enabled,
  so an intentionally incomplete intermediate camera is never presented as a
  one-frame jump while an operator watches the replay.
- Lazy hierarchy expansion, exact-path selection/highlighting, nested
  erase/redraw while selected, and selection clearing for models with a
  declared hierarchy probe.
- Production-main qged reference captures where an X11 display is available.
- Default or explicitly selected swap intervals.

The graphical scripts source `qged_test_display.sh`.  It accepts an existing
display only after `xdpyinfo` proves it is live and otherwise starts a private
Xvfb server.  This prevents a stale inherited `DISPLAY` from turning a product
test into an unrelated connection failure.  A checkpoint after a retained
scene mutation is an observation barrier: wait for progressive idle and an
actual presentation before sampling pixels.  OSMesa export performs that
presentation naturally, while a fast System-GL front buffer can otherwise
retain the preceding frame.

The model profiles are:

- `smoke`: Generic Twin and Lucy.
- `full`: smoke plus havoc and Hubble.
- `stress`: full plus shared Lucy, xpushed Lucy, and combined Stanford scenes.

Repeated-instance Lucy and distinct-asset scale are deliberately separate.
Generate the latter with:

```sh
ninja -C .build qged_unique_mesh_stress_generator
./.build/src/qged/qged_unique_mesh_stress_generator \
  .build/unique_mesh_stress.g
```

The default fixture has 5,000 independently stored, independently perturbed
BoTs spanning approximately 500 to 100,000 faces.  The current deterministic
distribution contains 4,962,125 source faces; its long tail reaches 100k
without making every leaf large.  Exactly 95 percent are closed, consistently oriented,
shared-index manifold meshes.  The remaining five percent is divided equally
between open surfaces and deliberately non-manifold closed meshes; leaf names
identify all three topology classes.  This ratio makes normal vehicle-part
behavior the primary correctness and performance signal while retaining
adverse-input coverage.

The deterministic shape profile is modeled after the visual distribution of
the Hubble qualification model rather than a uniform mesh array.  It combines
rounded main bodies, cylinders, long thin booms, box-like equipment, very thin
panels, open dishes, and irregular housings.  Physical extent varies
independently of triangle count across more than two orders of magnitude,
from hull-scale skins through structural pieces, subsystems, ordinary
components, and bolt-scale leaves.  Consecutive hierarchy branches form
recognizable bus, paired panel-wing, truss/boom, and instrument-cluster macro
assemblies.  This deliberately correlates hierarchy locality with different
visual significance: a prefix-biased realization or LoD allocator can no
longer appear representative merely because every branch contains the same
uniform sample.  Sixteen-leaf groups mix colored regions and ordinary
combinations, and balanced eight-way parent levels provide a deep hierarchy.
The generator rejects a 1,000-or-larger fixture if any shape family, the thin
aspect-ratio tail, or the required physical-size range is absent.
It therefore stresses source import, PoP construction, cache fan-out, result
publication, upload, hierarchy expansion, region/color handling, and memory
pressure rather than shared-asset instancing.  `--profile stress` includes it
when the database exists; it may also be selected explicitly as `--cases
unique_mesh_stress`.  Increase the count or maximum grid dimension for
release-candidate hardware characterization.

For an asset-count characterization an order of magnitude larger, generate a
separate 50,000-leaf database while capping each source mesh at a 24-cell grid:

```sh
./.build/src/qged/qged_unique_mesh_stress_generator \
  .build/unique_mesh_50k_stress.g 50000 24
src/qged/tests/qged_gui_matrix.sh --profile stress \
  --cases unique_mesh_50k_stress --backends system --modes shaded \
  --no-baseline --perf unique_mesh_50k_stress --timeout 900
```

Large cache characterizations should reuse that cache for follow-up profiling
instead of creating another multi-gigabyte copy.  `--warm-cache` deliberately
accepts only one case/backend/mode/swap combination and skips the cold phase.
It accepts only a cache marked by a previously validated cold matrix run for
the same canonical database file and file identity; merely pointing at a
nonempty cache directory is not proof that this model is warm:

```sh
src/qged/tests/qged_gui_matrix.sh --profile stress \
  --cases unique_mesh_50k_stress --backends system --modes shaded \
  --no-baseline --perf unique_mesh_50k_stress --perf-phase warm \
  --warm-cache /path/to/prior/caches/unique_mesh_50k_stress-system-shaded-swapdefault/cache \
  --timeout 900
```

This keeps the source-face population near 50 million so the comparison
isolates distinct-asset, hierarchy, scheduling, cache-file, and draw-submission
scaling rather than simply making the mesh construction workload ten times
larger.  It is explicit rather than part of the default stress profile because
its cold and warm caches are release-characterization artifacts.

For the realistic multi-gigabyte tier, generate 150,000 distinct assets at
the same bounded per-asset detail:

```sh
./.build/src/qged/qged_unique_mesh_stress_generator \
  .build/unique_mesh_150k_stress.g 150000 24
src/qged/tests/qged_gui_matrix.sh --profile stress \
  --cases unique_mesh_150k_stress --backends system --modes shaded \
  --no-baseline --perf unique_mesh_150k_stress --timeout 1800
```

The resulting database is expected to exceed 3 GB on disk and contains
142,500 closed manifold meshes, 3,750 open surfaces, and 3,750 deliberately
non-manifold meshes.  This tier is an explicit release/hardware qualification
case.  Its gates distinguish time-to-first-overall-extent, time-to-all-leaf
boxes, interactive latency during cold preparation, warm-cache mesh
publication, and eventual view-budget convergence; it must not be interpreted
as requiring all source triangles to be resident.

The 50k and 150k replays deliberately rotate the cold scene at the 1.5-second
first-useful checkpoint, before waiting for background work to settle.
First-useful pixels and cold input/render latency are therefore startup
requirements; eventual idle is a separate terminal characterization.  Camera
zoom/rotation then share one recovery wait instead of forcing complete cache
population after every synthetic input.  This distinction prevents a long
cache-preparation time from masquerading as unusable startup, while still
requiring the background state machine to terminate.

The smoke profile likewise gives terminal quiescence 30 seconds while keeping
its first-useful framebuffer gate and terminal non-box PoP-payload assertion
separate.  Lucy contains roughly 28 million source faces: on a genuinely empty
cache the final PoP cache write, bounded cut progression, calibration cooldown,
and resident compaction may finish after 15 seconds without implying that the
user waited that long for a useful visual.

Cold mesh realization is admitted by both worker count and an estimated
transient working set based on source face/point counts.  Tasks which have not
decoded those counts yet, including cold AABB and PoP-cache preparation, make a
conservative reservation from the serialized database-object size; a cheap
proxy result must not imply a zero-cost import.  The automatic byte allowance
is the lesser of 1 GiB, one-eighth of installed RAM, and one-quarter of
currently available RAM.  A task larger than the allowance runs alone;
`BObolLodService::setWorkingSetLimit` provides an explicit host/application
override.  Reports record the configured allowance, current and peak
reservations, current and peak executing-task counts, and retained mesh bytes
so parallel speedups cannot hide memory oversubscription.

The runner inventories each input database before use and verifies its size
afterward.  A suspiciously small `./.build/Generic_Twin.g` is rejected in favor
of the installed `share/db/faa/Generic_Twin.g`, preventing a blank database
from producing a misleading pass.  Each completed replay is also rejected if
its report records failed sources, payloads detached from compact entries,
superseded fallback presentations, missing retained draw roots, missing
selection state, or zero-length checkpoints.  The 1.5-second framebuffer is
also inspected for actual model-area pixels, so a retained but empty
database-source placeholder cannot satisfy the first-useful-frame gate.
The inspected model area excludes the border, bottom HUD, and right-side
convergence indicator; those widgets must never be counted as geometry.
Lucy's terminal check is stricter still: it requires a resident PoP asset,
nonzero progressive faces, an entry-backed CAD payload, filled model pixels,
and a satisfied prominent-object projected-error floor in the fixed shaded
`ae 90 0` view.  Cut ordinals are producer-defined and deliberately have no
quality meaning; accepting either an arbitrary minimum ordinal or merely a
nonzero mesh allowed premature compaction to masquerade as a successful
terminal view.  A correctly framed box followed by an empty, permanently
skipped, or unrecognizably coarse mesh is a failure in both cold and warm runs.
Its smooth-zoom probe keeps a gesture active while low-amplitude wheel events
arrive faster than the quiet debounce.  Before release, both backends must
load the missing pixel-demanded suffix and submit a richer effective PoP cut
than the pre-zoom stable image.  That proof may occur at the fixed
in-gesture checkpoint or the later low-amplitude checkpoint: OSMesa may use the
render-only ceiling to back down a subsequently measured discrete quality
probe at 100 ms, but the retained occurrence must not restart at its minimum
cut.  System GL additionally requires flat cumulative
full-upload bytes and increasing suffix-upload, lineage-reuse, and (when the
buffer grows) device-copy bytes.  The following quiet checkpoints wait past
the 750 ms memory-maintenance delay and prove that prefetched residency is
compacted without losing the useful same-scale return cut.

Quality checkpoints are passive samples.  The System-GL canvas reads its
already-presented Qt framebuffer without requesting another paint; the OSMesa
canvas uses an observational export traversal which does not advance
progressive admission or report a presentation-timing sample.  This is a
deliberate distinction from `QOpenGLWidget::grabFramebuffer()`, which can
request and time a render and thereby make capture frequency alter the LoD
cut.  Continuous-interaction checkpoints likewise observe the bounded state
while events are still arriving rather than forcing extra progress.

For hierarchy probes, framebuffer comparisons additionally require visible
whole-row tree selection, selected-subpath erase, redraw, and selection-clear
transitions.  Those hierarchy operations must return within 250 ms and may
not fall back to a broad tree scan.  Large OSMesa progressive scenes must
demonstrate a smaller face cut and a render below 250 ms while rotation is
held; their quiet terminal image remains a separate quality/performance
measurement.

`test_qged_measure_ui.json` is the focused measurement contract.  It uses a
closed, oriented BoT and drives qged's actual palette and canvas through 2D
and exact-hit 3D distance/angle gestures, degree/radian reporting,
right-click cancellation, resize, and replacement by another view tool.  The
single- and quad-layout CTests use OSMesa; the companion CMake runner accepts
`USE_SYSTEM_GL=ON` for an X11 qualification run and `ARTIFACT_DIR` to retain
its images.  The replay checks semantic overlay point counts and numeric
readback, then compares its own captures: a completed guide must alter the
framebuffer and cancellation must restore the exact baseline.  Measurement
line layers deliberately disable depth testing so a 2D guide through shaded
geometry remains visible; other line-layer overlays remain depth-tested by
default.

`test_qged_view_settings_ui.json` is the focused faceplate and view-settings
contract.  It drives the real palette and toolbar in both single and quad
layouts, then checks bidirectional GED/widget state for ADC, center dot, grid,
axes, scale, parameter/FPS text, and framebuffer composition.  The replay also
checks the corresponding retained Obol features, applies a world-space cutting
plane to shaded geometry, requires a visible image delta, and requires
disabling the plane to restore the exact baseline.  The companion CMake runner
uses OSMesa by default and accepts `USE_SYSTEM_GL=ON` and `ARTIFACT_DIR` like
the measurement runner.  Command-driven control refresh is keyed by libbv's
monotonic frame revision; a renderer-consumable dirty latch is not a sufficient
notification contract.

Typical runs:

```sh
src/qged/tests/qged_gui_matrix.sh --profile smoke --no-baseline

src/qged/tests/qged_gui_matrix.sh --profile full \
  --main-build-dir /home/cyapp/brlcad/.build

src/qged/tests/qged_gui_matrix.sh --profile stress \
  --perf stanford --apitrace generic_twin

src/qged/tests/qged_gui_matrix.sh --cases generic_twin \
  --backends system --modes shaded --swap-intervals default,0,1
```

qged's OSMesa backend currently uses Obol's retained fixed-function VBO path
by default.  It preserves the same scene/LoD policy and PoP cut selection as
the System GL renderer, but uses the faster software presentation kernel.
Both paths consume the same explicit ambient/material/light snapshot; the
software renderer does not inherit `GL_COLOR_MATERIAL` or stale fixed-function
light slots from an earlier frame.  Smoke runs for Generic Twin and Lucy also
capture the `mged` and `studio` lighting profiles and validate their ambient
intensity and camera-light counts.  Matching geometry, view, and PoP cuts are
expected to match visually across System GL and OSMesa; software execution is
allowed to be slower and its frame budget may select a coarser cut.
Set `OBOL_CAD_SOFTWARE_GLSL=1` to qualify the programmable OSMesa path through
the same cold/warm matrix.  The opt-in is intentional: the OSMesa GLSL
interpreter now has correctness-checked four-lane vertex and fragment
execution, but CAD scenes with many small triangles still spend substantial
time in programmable vertex dispatch and the generic span rasterizer.  Do not
make GLSL the software default until representative Generic Twin, Lucy, and
distinct-mesh stress runs meet or beat the fixed path without changing
framebuffer correctness.

The version-2 `report.json` beside each image set records render duration,
presentation
cadence/FPS, pending progressive work, LoD view/policy revisions, submitted and
applied work, active PoP cuts and face/point counts, cache-service
activity, realization failures, and diagnostics after every replay event.
Command text and checkpoint paths are recorded with their samples so
draw/erase phases can be compared without relying on brittle event indexes.
Full compact leaf invariants are collected only for the final sample unless
`QGED_TEST_DEEP_LOD_REPORT=1` is set, so transient sampling does not itself
distort large-model interaction timing.  This makes it possible to distinguish
“the screenshot eventually looks right” from an interactive system that
remained blocked, restarted PoP generation on view changes, or silently
stopped refining.
Set `QGED_TEST_DEEP_LOD_REPORT=0` for profiler runs to suppress even the final
full hierarchy inventory; at tens of thousands of occurrences that diagnostic
can otherwise contribute about a second of observer work unrelated to the
renderer being measured.

`--perf CASE` and `--apitrace CASE` apply only to the cold, shaded, default-swap
System OpenGL run of that case.  They are deliberately opt-in because tracing
changes timing and can generate large artifacts.
