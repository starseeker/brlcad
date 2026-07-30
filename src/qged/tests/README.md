# qged graphical validation

`qged_gui_matrix.sh` drives the real qged Qt canvas through deterministic event
streams.  It is intended for release-candidate and investigative validation,
not as a short CTest.  Each run writes screenshots, command logs, controller
timing/LoD samples, cache inventories, and a summary into a new artifact
directory.

The primary matrix dimensions are:

- System OpenGL and qged's `-s` OSMesa/software path.
- A newly created, empty `BU_DIR_CACHE` with no format or data children,
  followed by a warm process using the exact same cache.  (`BU_DIR_CACHE`
  itself must exist or libbu intentionally disables caching.)
- Shaded and wireframe PoP presentations.
- A deterministic 1100×800 top-level resize through the real Qt layout and GL
  resize path, independent of the operator's saved qged window settings.
- Initial draw return, 50 ms, 200 ms, 1.5 seconds, and stable checkpoints.
- `ae 90 0`, zoom in/out/return, sustained held-button rotation checkpoints,
  and post-release stable checkpoints.
- Lazy hierarchy expansion, exact-path selection/highlighting, nested
  erase/redraw while selected, and selection clearing for models with a
  declared hierarchy probe.
- Production-main qged reference captures where an X11 display is available.
- Default or explicitly selected swap intervals.

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
BoTs spanning approximately 500 to 100,000 faces, with an aggregate near 42
million source faces.  Exactly 95 percent are closed, consistently oriented,
shared-index manifold ellipsoids.  The remaining five percent is divided
equally between open heightfield surfaces and deliberately non-manifold closed
meshes; leaf names identify all three classes.  This ratio makes normal
vehicle-part behavior the primary correctness and performance signal while
retaining adverse-input coverage.  Physical extent varies independently of
triangle count and spans hull-scale skins through structural pieces,
subsystems, ordinary components, and bolt-scale leaves.  The leaves occupy a
three-dimensional elongated envelope with varied orientations.  Sixteen-leaf
groups mix colored regions and ordinary combinations, and balanced eight-way
parent levels provide a deep hierarchy.
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
cache the final PoP cache write, bounded level walking, calibration cooldown,
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
nonzero progressive faces, an entry-backed CAD payload, and filled model
pixels.  A correctly framed box followed by an empty or permanently skipped
mesh is a failure in both cold and warm runs.
For hierarchy probes, framebuffer comparisons additionally require visible
whole-row tree selection, selected-subpath erase, redraw, and selection-clear
transitions.  Those hierarchy operations must return within 250 ms and may
not fall back to a broad tree scan.  Large OSMesa progressive scenes must
demonstrate a smaller face cut and a render below 250 ms while rotation is
held; their quiet terminal image remains a separate quality/performance
measurement.

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
Set `OBOL_CAD_SOFTWARE_GLSL=1` to qualify the programmable OSMesa path through
the same cold/warm matrix.  The opt-in is intentional: the OSMesa GLSL
interpreter now has correctness-checked four-lane vertex and fragment
execution, but CAD scenes with many small triangles still spend substantial
time in programmable vertex dispatch and the generic span rasterizer.  Do not
make GLSL the software default until representative Generic Twin, Lucy, and
distinct-mesh stress runs meet or beat the fixed path without changing
framebuffer correctness.

The `report.json` beside each image set records render duration, presentation
cadence/FPS, pending progressive work, LoD view/policy revisions, submitted and
applied work, active PoP cut levels and face/point counts, cache-service
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
