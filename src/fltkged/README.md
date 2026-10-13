# fltkged

`fltkged` is the FLTK MGED application built on `libfltkcad`. It supplies the
menu and application lifecycle around the library's terminal, hierarchy,
editing, and display widgets, allowing work with a BRL-CAD database without a
Tk dependency. When Tcl is available
at build time, a bounded embedded interpreter supports scripts and dispatches
GED commands. Legacy Tk event servicing is available on demand, while new
scripts can open native panels through the `fltkged` Tcl command.

The default layout uses separate command and graphics windows. Start with
`--combined` or select Settings > Combined Windows to place graphics above the
terminal in one window. Tools > Database and View Controls opens the hierarchy,
attributes, camera, and interaction panels. Menus are available in both main
windows; closing either exits the application after checking unsaved edits.

The initial graphics view uses MGED's azimuth 35, elevation 25 orientation.
Drawing into an empty display automatically frames the geometry; subsequent
draws preserve the camera. `draw -R` suppresses automatic framing, and `B`
replaces the display and frames its new contents.

`clear` clears the terminal display while retaining command history and
geometry. `quit`, `q`, and `exit` use the same unsaved-edit checks as the File
menu. Database commands (`opendb`/`open` and `closedb`/`close`) update the native
panels and refuse to discard unsaved edits. `search -exec` dispatches GED
commands and sends their output to the terminal. These behaviors also apply
to GED commands dispatched through the embedded Tcl bridge.

View size, center, and AET are libdm overlays inside each graphics view.
Modes > Faceplate (F7) toggles them. Modes > Faceplate GUI (F8) enables the
optional in-view button menu and sliders. Middle-click a menu entry or drag
a slider to adjust the view. These controls use the same toolkit-independent
libdm renderer and hit testing as MGED; the application supplies their actions.
The FLTK menu supplies camera, stored-view, and ADC controls while viewing.
During editing it displays the shared primitive or instance menu together with
Accept, Apply, Reset, and Reject. Both frontends supply these actions through
the libdm callback interface.

The graphics window uses MGED's default mouse bindings:

| Gesture | Action |
| --- | --- |
| Left / right click | Zoom out / in by two |
| Middle click or drag | Pick, query, or sweep using the selected Modes > Mouse tool |
| Shift + any button drag | Translate the view |
| Control + any button drag | Rotate the view |
| Shift-Control + any button drag | Scale the view |
| Alt with those modifiers | Constrain to X / Y / Z using left / middle / right |
| Mouse wheel | MGED wheel zoom factors and direction |

Help opens the local HTML manual browser, with section and name filters, page
search, and Back/Forward navigation. In FLTK builds, `brlman --gui` uses this
same viewer. `brlman --gui mged` opens a named page; `brlman --gui -S n draw`
selects a GED command page. Enable `BRLCAD_EXTRADOCS` to build and install the
HTML manuals.

## Implemented interface

- single and four-view display layouts using either system OpenGL or OSMesa;
- a bounded command terminal with GED continuation prompts, history,
  completion, and asynchronous subprocess output;
- a lazy database hierarchy with draw/erase actions, object attributes, and a
  bounded name/type search dialog;
- active-view camera and faceplate controls, including grid/snap, adaptive
  plotting, axes, framebuffer mode, and parameter overlays;
- primitive creation using libRT's current make-type registry;
- descriptor-driven primitive editing with in-memory preview, undo, and
  explicit save;
- transactional combination editing for the Boolean tree, region fields,
  color, shader, material, and inheritance;
- a dedicated sketch editor with direct vertex dragging, geometric constraints,
  line/arc/Bezier/NURB topology and control editing, preview, undo, and save;
- point/rectangle selection, query-ray, rectangle zoom, and object transforms;
- object and combination-member matrix editing;
- embedded and file raytracing, framebuffer IPC and image transfer tools;
- validated persistent display and terminal preferences;
- database creation, import, and selected-object export;
- validated BOT/NMG repair and material/density import, export, validation,
  source, list, set, and clear workflows;
- embedded HTML documentation plus an about window; and
- view presets, per-viewport stored-view rings, autoview, redraw, and database
  replacement.

Reusable widget headers are installed under `fltkcad/`; `fltkged` retains only
its application window and localization code. Editors keep changes in memory
until Save. They detect deleted or externally replaced database objects before
writing, and database replacement or application exit cannot silently leave a
live edit session referring to a closed database.

Classic editing commands and native primitive/transform panels use the same
libged session. `sed path` edits a primitive definition and `oed path` edits an
instance matrix. `accept` saves and ends the edit; `sed_apply` saves and keeps
editing; `sed_reset` reloads the saved geometry; `reject` discards changes.
Working geometry appears in the main viewports before it is saved. Transform
drags and panel actions retain this working state until Accept or Reject.
Primitive operations use librt's menus and parameter descriptors. Model/view/
object coordinates, keypoints, parameter and normalized mouse input, and knob
rates are available through the shared command interface. FLTK advances knob
rates from its event loop. Tk widget paths and synthetic Tk mouse events remain
native MGED interfaces.

On POSIX systems, GED subprocess streams and framebuffer clients are integrated
with `Fl::add_fd`. Windows subprocess pipes use bounded FLTK timer polling and
BRL-CAD's nonblocking process API. Output is appended by delta and bounded
before it reaches the terminal. Active subprocesses, watches, and framebuffer
clients are cancelled during database replacement or window teardown.

## Failure and performance behavior

The application reports missing databases, unsupported object types,
read-only databases, stale edit targets, invalid numeric values, backend
failures, and subprocess errors through the relevant widget and terminal.
Unsaved editor changes must be saved or deliberately discarded before opening
another database or leaving the application.

Potentially unbounded inputs are capped: command input and accumulated
arguments, GED output, completion lists, hierarchy roots and children, search
results, combination text/tree data, and sketch vertex/segment counts all have
explicit limits. Hierarchy children load only when a combination is expanded.
Primitive and sketch previews operate on transient in-memory internals and do
not write the database during interactive changes.

Tcl scripts and results are limited to 4 MiB, each evaluation has time and
command-count limits, and the Tcl event pump handles at most 32 queued events per
FLTK timer tick. Density files are limited to 64 MiB; imports are validated
first and require confirmation because BRL-CAD replaces the current density
table. Density export and geometry-copy operations refuse existing outputs.

The HTML viewer and terminal apply their own document/resource, decoded-image,
history, and scrollback limits; see `../libfltkcad/docs` for exact defaults.

## Build and test

`fltkged` is built as part of BRL-CAD. From the BRL-CAD source root:

```sh
cmake -S . -B build \
  -DBRLCAD_EXT_DIR=/path/to/bext/build \
  -DBRLCAD_ENABLE_FLTK=ON
cmake --build build
ctest --test-dir build -R '^(fltkcad|fltkged)_' --output-on-failure
```

FLTK and litehtml must be installed by bext. `BRLCAD_COMPONENTS=fltkged`
selects the editor and its dependencies. The executable installs into BRL-CAD's
binary directory alongside its other applications. Tcl scripting follows
`BRLCAD_ENABLE_TCL`; the interface can also be built with Tcl disabled.

The command prompt defaults to Tcl when Tcl is enabled. Substitutions,
procedures, loops, `source`, multiline commands, and `puts` work directly at
the prompt. Incomplete Tcl input and GED commands requesting more arguments
retain a continuation prompt. Settings / Command Language switches between
Tcl and raw libged arguments; `--shell=tcl` and `--shell=ged` select the startup
language. Switching languages preserves the Tcl interpreter's state. A build
without Tcl uses the raw prompt.

`fltkged -c database.g 'source {script.tcl}'` runs without a display server;
`fltkged -c database.g < script.tcl` reads Tcl commands from standard input.
Batch mode creates a missing database, accepts MGED's `opendb file.g y`
creation form, and supports `quit`, `q`, and Tcl `exit ?status?`.
`--shell=ged` also works in batch mode. Native menus and editor panels continue
to issue GED commands directly, independently of the prompt language. The
terminal's command-handler callback is the boundary for other interpreters;
additional languages are not bundled.

Tcl timers and file events run while FLTK owns the main loop. A Tcl event
source also services FLTK while scripts wait in `vwait` or `update`, without
replacing Tcl's process-wide notifier. Event servicing uses a bounded 10 ms
poll. Tk remains opt-in. Scripts retain the existing input, result, time, and
command-count limits; the command budget is renewed for each evaluation.

See [`REGRESSION_STATUS.md`](REGRESSION_STATUS.md) for measured MGED test
compatibility, coverage distinctions, and the isolated audit runner.

Headless tests run under Xmin and cover startup, missing databases, both view
backends, embedded framebuffer raytracing, primitive/combination/sketch
editing, search, creation, Tcl dispatch, transforms, interaction modes, and
active-view controls. The source database is copied to a temporary file before
any editor fixture is created.

See [`COMPATIBILITY.md`](COMPATIBILITY.md) for the explicit MGED parity target,
implemented native workflows, and compatibility policy. Human release checks
are defined in [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md).

Preferences includes terminal foreground/background RGB colors and a status
font field accepting `auto` or a fixed size. Defaults are black terminal text
on white and automatic status sizing. Apply and Save persists these settings.
The status text stays in the graphics viewport and wraps as needed.
