# MGED regression compatibility

The FLTK/Tcl/Tk-enabled configuration passes **119 of 119 tests that
actually invoke MGED**, with `fltkged` substituted for every MGED invocation.
The same isolated audit using native MGED passes **119 of 119**.

| Coverage | fltkged | MGED |
| --- | ---: | ---: |
| MGED command scripts | 85 / 85 | 85 / 85 |
| Regression harness failure-detection checks | 1 / 1 | 1 / 1 |
| Other integration and regression tests using MGED | 33 / 33 | 33 / 33 |
| Total | 119 / 119 | 119 / 119 |

Twelve additional IGES test configurations accept an MGED executable but do
not invoke it for their selected cases. They are excluded from these totals.
The two fixture cleanup tests are also excluded because they run CMake only.
Counts depend on the configured components and describe the assertions that
ran, rather than complete command or GUI equivalence.

Passing coverage includes database construction, combination operations,
classic primitive and instance editing, view commands, Tcl helpers, sourced
procedures and aliases in `search -exec`, asynchronous overlap plotting, NIRT,
BOT operations, DSP, IGES, NACA geometry, red, tree, weight, solid rendering,
and the shader image comparison. These tests exercise the real FLTK
executable's headless command path; they do not launch MGED as a helper.

## Shared editing implementation

The previous comparison was 92/119 for FLTK and 118/119 for MGED. The remaining
editing failures are now addressed through a per-GED session in libged, with
its public API in `ged/edit.h`. MGED borrows its persistent librt edit context;
FLTK's prompt, primitive panel, transform panel, object drags, and faceplate
menus use the same session. Toolkit widgets and Tcl interpreters stay outside
that session.

The shared state handles selection and repeated-member identity, primitive
working geometry, instance matrices, keypoints, coordinates, contextual
transforms, parameter and normalized mouse input, knob rates, and
accept/apply/reset/reject. Preview vectors are temporary; the database changes
only on save. Failed saves preserve working geometry, and external changes or
deletion of a selected path are detected before writing. Reset retains host
callbacks. Scripted `edit -i` operations on the active primitive use the same
working geometry. Database closure clears both interactive and standalone
scripted buffers before releasing the database.

Primitive operations use librt handlers, menus, and descriptors. ARB face
definition is shared by MGED and libged through `rt_arb_facedef`, replacing the
duplicated implementations. Common ADC commands, standard views, status unit
fields, and Tcl unit-variable traces cover the other former state failures.
Native MGED now recognizes the documented `status base2local` and
`status local2base` names as well as its old internal aliases.

The shader fixture explicitly establishes its initial projection camera,
removing its dependence on MGED's startup and implicit drawing behavior.
The full shader pixel comparison passes in both frontends.

## Coverage boundaries

`csg_edit` runs shared geometry and edit-command assertions in both frontends.
Native MGED also runs its selection and synthetic display-manager mouse
blocks. FLTK reports their exclusion with a `COVERAGE` message and exercises
the common edit operations directly. Tk widget paths, `_mged_matpick`, and
synthetic X events are not portable interfaces. The tests retain native
coverage rather than pretending those events ran in FLTK.

The command inventory test checks FLTK's installed documented GED commands
and help, using `::ged::commands`. Native MGED retains its application command
inventory and no-argument invocation checks. These modes have different
coverage: ordinary libged commands report missing-argument errors, and native
MGED has application-specific commands. Both modes now wait for the final
partial batch of parallel checks, which was previously left unchecked. Tcl
help is loaded before GED registration so its `help` and `?` procedures are
preserved; raw GED help remains accessible in the GED namespace.

The shared `rset` interface currently covers coordinates and view/edit
transform selection. Shared `status` exposes state, scale, units, and view/edit
matrices; its no-argument result is the state rather than MGED's full display
summary. ADC-specific GUI gestures and other native application resource
settings remain separate work. Default FLTK modifier drags manipulate the
view; explicit object tools edit working instance matrices.

Additional tests verify session isolation, read-only and stale saves, callback
retention, switching active views, scripted working-copy reuse, knob ticks,
repeated members under parent transforms, nonuniform scaling, transient
viewport previews on both rendering backends, Tcl unit rebinding, edit menus,
and database/quit guards. A passing script remains only as precise as its
assertions; some older scripts are still smoke tests.

## Validation

The full build and 1,260 configured CTest tests pass, including the shared
editing unit test, FLTK application guards, both rendering backends, Tcl/Tk
integration, and native MGED regressions. This run excludes the informational
`NOTE:` entry, the expected-failure `bu_heap_1` test, and `slow-bigdb_5gb`;
the 1 GB database stress test is included. The isolated FLTK audit passes all
119 tests above. `regress-repository` and `regress-licenses` also pass.

## Corrected tests and stronger checks

The six previously failing scripts now pass in both frontends. Their fixes
correct clone naming expectations, use direct parent/child paths for `rm`,
supply the missing `killtree` union operator, avoid reading deleted objects,
close the opening `hide` string, and prepare the saved-view and plot inputs.

All command scripts are now sourced through a driver that catches Tcl errors,
requires its own completion marker, and rejects explicit assertion failures.
The harness checks that errors, incomplete input, failed assertions, printed
failures, and premature quit/exit cannot pass, even when the frontend returns
success. Geometry patterns are expanded explicitly rather than depending on
MGED's interactive glob pass, which had escaped Tcl substitutions in assertions.

New assertions check clone translations, hierarchy references, pattern copies,
safe member removal, deletion, visibility, primitive creation, restored views,
and overlay display vectors. Orthographic view restoration checks projected
center rather than requiring an unchanged depth coordinate: saveview backs
the ray origin out of the geometry along the view axis.

Each test has its own working directory. CTest fixtures order `saveview` before
`loadview` and `plot` before `overlay`, with cleanup after the consumers. The
audit honors the same setup, failure propagation, and cleanup ordering. The
export tests now use `postscript`; MGED's `ps` command lists processes.
The `rotobj` test uses supported absolute/incremental rotations, and the
search output-hook test emits actual `bu_log` output inside sourced scripts.

## Reproduce the comparison

Build BRL-CAD with FLTK and Tcl enabled, then run from the source tree:

```sh
python3 src/fltkged/tests/regression_audit.py ../build --output-directory ../audit-fltkged
python3 src/fltkged/tests/regression_audit.py ../build --frontend mged --output-directory ../audit-mged
```

The audit discovers configured CTest commands, substitutes the frontend,
provides a private `ensearch mged` wrapper, and records actual invocations.
Each run has a fresh workspace with its own databases, caches, logs, and
JSON results. Selecting a fixture consumer includes its setup and cleanup.
The installed build executables are not replaced. Tests that
never invoke the wrapper are reported separately. The audit returns nonzero
when any test fails.
Use `--include '^regress-mged_'` to select the command suite, or another test
name expression to investigate a smaller group.

The audit currently requires a POSIX host and an out-of-source build without
an alternate `mged` in the source tree's `bin` directory. Additional enabled
components may contribute more tests than the configuration measured here.
