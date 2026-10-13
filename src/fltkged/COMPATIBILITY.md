# MGED compatibility matrix

This matrix defines the intended replacement boundary for `fltkged`.  A
feature marked **GUI** has a native FLTK workflow, **Command** is available
through the terminal or embedded Tcl, **Compatibility** denotes an optional
legacy adapter, and **Validation** denotes a release gate rather than missing
application code.

| MGED capability | Status | FLTK implementation and limits |
| --- | --- | --- |
| Open, create, import, and export databases | GUI | Existing files open read/write; creation and export refuse accidental overwrite; import uses `concat` conflict handling. |
| Window arrangement | GUI | Separate command/graphics windows by default, optional combined panes, and a separate database/view controls window. |
| Command window | GUI | Editing in the terminal display, UTF-8 cursor operations, clipboard input, interrupt handling, bounded scrollback, command history, completion, continuation prompts, and asynchronous GED process output. |
| Tcl scripts | GUI/Command | Default prompt language when enabled; substitutions, procedures, multiline input, and sourced scripts share `ged`/`::ged::*` dispatch and the native-panel bridge. Raw GED is selectable. Tcl/FLTK events cooperate without requiring Tk. |
| Batch execution | Command | `-c` accepts command arguments or Tcl input on stdin without a display server. See `REGRESSION_STATUS.md` for the measured compatibility boundary and exit-status differences. |
| Geometry display | GUI | Single or four views, system OpenGL or OSMesa, view presets, mouse navigation, capture, and bounded framebuffer sizes. |
| Embedded framebuffer | GUI | `fbserv` IPC, progressive `ert`, clear/import/export tools, overlay and underlay modes on both display backends. |
| Picking and interaction modes | GUI | MGED click zoom and modifier drags on all three buttons, axis constraints, middle-button tools, point selection, rectangle selection/zoom, query ray, and object translate/rotate/scale. |
| View and faceplate controls | GUI | AET, center, size, perspective, grid and snapping, ADC, axes, scale, z clipping, adaptive plotting, parameters, and framebuffer mode. |
| Stored views | GUI | A ten-entry ring per viewport for the current session. Named database view objects remain available as GED commands. |
| Object hierarchy and attributes | GUI | Lazy bounded hierarchy, draw/erase, selection, attribute editing, and bounded database search. |
| Primitive creation | GUI | Uses BRL-CAD's current make-type registry rather than a duplicated primitive list. |
| Primitive editing | GUI | Descriptor-driven `libged`/`librt` operations with preview, undo, stale-object detection, and explicit save. Classic type-specific forms are not duplicated. |
| Object and member matrices | GUI/Command | Shared libged working matrices, occurrence-aware member selection, keypoints, model/view/object coordinates, uniform and axis scaling, preview, and accept/apply/reset/reject. |
| Combination editing | GUI | Boolean tree, union/subtract/intersect insertion, region fields, color, shader, material, inheritance, reload, and explicit save. |
| Sketch editing | GUI | Direct pointer manipulation, vertex-list motion, coincident/horizontal/vertical/distance constraints, tangency, splitting/reversal, arc and NURB controls, preview, undo, and explicit save. |
| Raytrace and query ray | GUI | Embedded render, file render, abort, center query, framebuffer clear, `png2fb`, and `fb2pix`. |
| Analysis tools | GUI/Command | Selected-object launchers for `analyze`, `gqa`, and overlap checking; their full option sets remain available in the terminal. |
| Help and about | GUI | FLTK `brlman` and a shared local HTML viewer with a filtered manual index, navigation history, page search, and missing/invalid-document diagnostics and a native about window. |
| Preferences | GUI | Persisted display gradient and terminal history/scrollback limits, with validation of corrupt or extreme stored values. |
| General GED command set | Command | Commands are dispatched directly through `libged`; multiline `GED_MORE` handling and bounded output are supported. |
| Legacy MGED Tcl/Tk GUI scripts | GUI/Compatibility | Tk is opt-in rather than a startup dependency. A bounded Tcl-event pump services legacy Tk windows, while `fltkged show`, `selection`, and `view` expose native panels to new scripts. Scripts coupled to a particular historical MGED widget path may still require adaptation. |
| NMG/BOT repair and material-library tools | GUI | Type-checked BOT inspection, repair, decimation, condensation, fusing and normals; NMG normals, collapse, simplify and labels; material and bounded density import/export/validation/edit operations. Destructive density changes require confirmation, and output-producing operations reject collisions. |
| Windows asynchronous subprocess integration | GUI | POSIX uses descriptor callbacks; Windows uses bounded FLTK timer polling through `bu_process_pending` and `bu_process_poll`, preserving the same libged callback lifecycle. |
| Accessibility, localization, and production UX review | Validation | Visible focus and control-name audits are automated, primary menus use gettext, and the required human platform matrix is defined in `RELEASE_CHECKLIST.md`. Human sign-off remains a release gate. |

## Compatibility policy

The terminal and Tcl bridge are deliberate compatibility layers. Legacy Tk is
loaded only when explicitly requested or required by a sourced script, so the
normal application remains an FLTK interface. A native panel is warranted when
it adds interactive display behavior, structured validation, transactional
editing, or safer file handling. Commands that are already clear and complete
at the GED prompt may remain command-driven.

Inputs that can grow with database or user data are bounded. Editors keep
changes in memory until an explicit save, verify that their database object has
not been replaced, and reject unsupported or stale targets. Closing a database
first terminates subprocesses and framebuffer clients, detaches views, and
then releases the GED context.

The classic solid/object editing lifecycle and commands share libged state
with native panels. Tk widget paths, native MGED selection automation, and
ADC-specific mouse gestures remain toolkit-specific. Current object transforms
use the explicit interaction tools and editor panels; the default modifier
bindings manipulate the view. A passing command regression does not establish
equivalence of every native GUI gesture.
