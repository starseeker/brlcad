# Main integration, September 2026

Status: rebased baseline with focused build and behavioral qualification.
Annotation image parity and the broader production gates remain open.  The
results below distinguish command behavior, retained rendering, GUI checks,
and evidence that still needs investigation.

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
- State-publication and an isolated full source-publication sweep pass.  An
  earlier source sweep failed the scene-light input assertion during concurrent
  qualification.  The isolated selector and subsequent full sequence pass;
  the intermittent result is retained under CONTROL-01, not declared repaired.
- Repository and license checks pass in the fresh original-checkout build.
  Loader checks resolve the recorded native dependency pair without an
  environment override. Reconfiguration retains that selection.

## Remaining limits and next acceptance

All seven annotation image comparisons differ from upstream controls.  The
confirmed duplicate-translation defect is repaired, but framing and rendering
still differ.  No image controls or comparison thresholds were replaced.
Qualify screen-space placement, text/fill/style rendering, bounds and autoview
before claiming full annotation visual parity (QUALITY-01 and EDIT-01).

The GUI checks use Xvfb and software rendering; they do not qualify native GPU,
Windows, macOS, fractional DPR, or every client workflow.  Physical cancellation
and close during active asynchronous commands, view unsharing, and the wider
large-model matrix remain required.  Passing the sourced search-exec regression
is not evidence for all interactive interruption cases.

The integration uses a recorded existing native Obol/OSMesa pair and preserves
its dirty-source patch and binary hashes.  This is not a clean rebuild or full
qualification of the native dependency sources.  BUILD-01 remains open.
The original checkout's new build directory is `.build-main`, with an isolated
prefix in `.build-main-deps`; the older `.build` is preserved and must not be
used as evidence for the rebased sources.  Qualification logs, images, conflict
ledger and dependency manifest are retained under
`.build-main/obol-qualification/20260919-main-integration`.

The broader production finish line remains
[the simplification guide](obol_simplification_guide.md).  Passing this
integration's checks establishes a usable new baseline; it does not close the
remaining resource, visual-quality, scale, or platform gates by itself.
