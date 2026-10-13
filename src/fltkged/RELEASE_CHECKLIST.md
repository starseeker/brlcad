# fltkged release validation

Automated tests are necessary but cannot replace human desktop review. A release
candidate is ready only after every applicable item below has an owner, a date,
and a recorded result. Failures belong in the release issue rather than being
silently waived.

## Supported platforms

- Exercise a release build on each supported Windows, Linux, and macOS version.
- Test system OpenGL and OSMesa where each backend is supported.
- Start and abort an embedded raytrace, then close and reopen a database while
  output is active. On Windows, confirm stdout and stderr remain responsive.
- Exercise framebuffer import, export, overlay, underlay, and resize.
- Confirm installation into a clean prefix and launch without build-tree paths.

## Keyboard and assistive use

- Reach every menu, tab, input, list, button, view, and terminal using only the
  keyboard. Confirm focus is always visible and follows a useful order.
- Verify menu mnemonics and shortcuts do not conflict on the target desktop.
- With the platform screen reader enabled, verify that labeled controls and
  status/error messages are announced intelligibly.
- At 100%, 150%, and 200% display scaling, confirm that text is not clipped and
  dialogs remain operable at the minimum supported screen size.
- Check focus, selection, wireframe, framebuffer, and error-state contrast in
  both the default palette and a high-contrast desktop theme.

## Localization

- Build and install at least one gettext catalog for the `fltkged` domain.
- Run under a non-English UTF-8 locale and a right-to-left locale.
- Confirm translated menu paths retain their slash hierarchy and mnemonic
  markers, long labels resize acceptably, and file/database names round-trip.
- Confirm terminal and Tcl input remain commands rather than translated prose.

## Modeling workflows

- Perform create/edit/save/reopen cycles for representative primitive types,
  combinations, member matrices, and sketches.
- Exercise direct sketch vertex dragging, each vertex constraint, tangency,
  splitting/reversal, arc changes, and NURB knot/weight changes.
- Inspect and repair representative valid and damaged BOT and NMG objects.
- Assign shaders/colors, edit density records, and verify undo or output-copy
  behavior before destructive specialist operations.
- Run representative non-GUI Tcl scripts. If legacy Tk support is shipped,
  run the supported legacy script set and close all Tk windows before exit.

## Failure and scale

- Repeat the documented missing, unreadable, malformed, stale, and read-only
  input cases for databases, HTML, Tcl, and framebuffer files.
- Exercise near-limit HTML documents, terminal history, GED output, hierarchy
  breadth, search results, combination trees, and sketches while monitoring
  responsiveness and memory growth.
- Confirm every rejected operation leaves database geometry unchanged and
  reports a useful error in the originating panel or terminal.
