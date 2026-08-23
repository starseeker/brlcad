# libBObol platform and threading matrix

Last reviewed: 2026-08-23

This is the supported-host statement for the direct libged/libBObol drawing
path.  “Deferred” means the core API remains portable, but the native host is
not represented as release-validated.

| Configuration | Status | Validated contract | Enabling work if deferred |
|---|---|---|---|
| Linux, no `DISPLAY` | Supported | Core scene/source, value-handle, LoD worker, retained RT, OSMesa, diagnostic/none endpoint, and headless-host tests | — |
| Linux/X11, Qt GL and software hosts | Supported architecture; current System GL matrix pending | Endpoint create/replace/detach, physical resize, input routing, framebuffer rehost, HW/SW/RT switching, and teardown | Re-run exact-current graphical matrix on a usable X host |
| Linux/X11, Tk GL and photo hosts | Supported architecture; current qualification pending | MGED/TclCAD endpoint lifecycle, input routing, framebuffer composition, pane deletion, and teardown | Re-run exact-current MGED/TclCAD matrix on a usable X host |
| Windows core/headless | Deferred | Source builds retain platform-neutral contracts | Validate OSMesa/headless package, thread and sanitizer jobs on MSVC/clang-cl |
| Windows WGL/Tk or Qt host | Deferred | WGL bridge is compiled where enabled | Add native CI host, resize/input/capture lifecycle tests, and validate endpoint replacement |
| macOS core/headless | Deferred | No release-support claim | Add supported Coin/Obol headless context manager and sanitizer CI |
| macOS native/Qt/Tk host | Deferred | No release-support claim | Select a native provider, implement factory registration, and add full lifecycle/capture tests |

## Owner-thread boundary

Each scene controller, view controller, display endpoint, native host, feature
store, polygon store, database source, and Coin scene graph has one owner
thread.  Mutation and queries that traverse or return borrowed scene state run
on that thread.  Endpoint callbacks may enqueue work or request a future frame;
they may not destroy their active controller reentrantly.

LoD, retained-RT, and detached database workers receive copied request values,
immutable geometry payloads, or database snapshots.  They return plain result
data.  They do not create, mutate, traverse, ref/unref, or publish Coin nodes,
touch a native window/OpenGL context, or initialize global Inventor types.  The
owner thread drains results, validates request/source revisions, and publishes
them.  Focused LoD and RT tests compare worker storage/thread identities with
the owner-thread publication identity.

## Hosted-view rule

Every intentional pane—including a hidden quad pane—retains one GED view
context, one display endpoint, and one controller.  No hidden pane shares a
controller or native host.  If lazy pane creation is introduced, pane
destroy/recreate coverage is required before it becomes a supported lifecycle.

The documented headless validation command is:

```sh
env -u DISPLAY -u WAYLAND_DISPLAY ctest --test-dir .build \
  --output-on-failure -L bobol_headless
```

The focused sanitizer label covers endpoint teardown, repeated attach/detach,
stale GED and feature handles, compact edit promotion/demotion, retained-RT
publication, and LoD worker shutdown.  Use separate build trees because TSan
cannot be combined with ASan:

```sh
cmake -S . -B .build-asan -G Ninja \
  -DBRLCAD_OPTIMIZED=OFF -DBRLCAD_DEBUGGING=ON \
  -DBRLCAD_ENABLE_ADDRESS_SANITIZER=ON \
  -DBRLCAD_ENABLE_UNDEFINED_SANITIZER=ON
ninja -C .build-asan
env -u DISPLAY -u WAYLAND_DISPLAY ctest --test-dir .build-asan \
  --output-on-failure -L bobol_sanitizer

cmake -S . -B .build-tsan -G Ninja \
  -DBRLCAD_OPTIMIZED=OFF -DBRLCAD_DEBUGGING=ON \
  -DBRLCAD_SANITIZE_THREAD=ON
ninja -C .build-tsan
env -u DISPLAY -u WAYLAND_DISPLAY ctest --test-dir .build-tsan \
  --output-on-failure -L bobol_sanitizer
```

ThreadSanitizer also depends on a host address-space layout that its runtime
can reserve.  In restricted or ptrace-based containers it may terminate before
`main` with `FATAL: ThreadSanitizer: unexpected memory mapping`.  That is a
runner limitation, not a test skip or a successful race check: such a runner
must still build the complete `bobol_sanitizer` target set, while runtime TSan
validation is performed on a compatible native Linux worker.  LeakSanitizer
has the analogous ptrace restriction; disabling leak detection is acceptable
only for ASan/UBSan diagnostics in that environment, not as the release leak
gate.
