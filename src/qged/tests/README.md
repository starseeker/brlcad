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

The runner inventories each input database before use and verifies its size
afterward.  A suspiciously small `./.build/Generic_Twin.g` is rejected in favor
of the installed `share/db/faa/Generic_Twin.g`, preventing a blank database
from producing a misleading pass.  Each completed replay is also rejected if
its report records failed sources, payloads detached from compact entries,
superseded fallback presentations, missing retained draw roots, missing
selection state, or zero-length checkpoints.  The 1.5-second framebuffer is
also inspected for actual model-area pixels, so a retained but empty
database-source placeholder cannot satisfy the first-useful-frame gate.
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

`--perf CASE` and `--apitrace CASE` apply only to the cold, shaded, default-swap
System OpenGL run of that case.  They are deliberately opt-in because tracing
changes timing and can generate large artifacts.
