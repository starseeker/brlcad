# libfltkcad

`libfltkcad` is the reusable FLTK interface library for BRL-CAD. It follows
the same application/library boundary as `libqtcad`: applications own their
top-level windows and policy, while this library supplies CAD-aware widgets.

The installed CMake target is `BRLCAD::fltkcad`, and public headers are
installed below `include/brlcad/fltkcad/`. The library includes:

- the `Fl_Cad_Terminal` command console and bounded command history;
- the `Fl_Cad_View` OpenGL/OSMesa display-manager widget and quad-view host;
- hierarchy, attribute, interaction, view-control, render, and search widgets;
- transactional primitive, combination, sketch, transform, BOT/NMG, and
  material editors;
- bounded Tcl/GED integration with optional legacy Tk event servicing; and
- the `Fl_Litehtml_View` local documentation viewer.

`fltkged` is the reference application using these widgets. It is deliberately
outside the library and owns the main window, menus, localization domain,
application preferences, and process lifecycle.

Primitive and transform widgets use the per-GED session in `ged/edit.h`.
`Fl_Cad_View::ged` binds a borrowed GED context for transient working-geometry
previews; clear that binding before closing the context. Transform panels and
object drags retain edits until Accept or Reject. The Tcl bridge dispatches
contextual editing commands through the same session and exposes current
database unit conversions without copying application editing state.

## Build and consume

Configure from the BRL-CAD source root with FLTK 1.4 or newer and litehtml
installed by bext. FLTK must include its OpenGL and image libraries and SVG
support. Static dependencies must be built with
`CMAKE_POSITION_INDEPENDENT_CODE=ON` for the shared library.
Dependency discovery uses BRL-CAD's standard bext lookup, preferring the
staged copies in the BRL-CAD build directory.

```sh
cmake -S . -B build \
  -DBRLCAD_EXT_DIR=/path/to/bext/build \
  -DBRLCAD_ENABLE_FLTK=ON
cmake --build build
cmake --install build --prefix /path/to/prefix
```

An installed consumer uses:

```cmake
find_package(BRLCAD CONFIG REQUIRED COMPONENTS fltkcad)
target_link_libraries(my_app PRIVATE BRLCAD::fltkcad)
```

```cpp
#include <fltkcad/Fl_Cad_Terminal.H>
#include <fltkcad/Fl_Cad_View.H>
```

`BRLCAD_ENABLE_FLTK` builds the library and `fltkged`; it defaults to `OFF`.
`BRLCAD_ENABLE_MINIMAL` disables FLTK support. For a component build, select
`BRLCAD_COMPONENTS=fltkged` or `BRLCAD_COMPONENTS=libfltkcad`.
`FLTKCAD_BUILD_EXAMPLES` defaults to `OFF`; `BUILD_TESTING` controls tests.
The test database defaults to `src/libged/tests/draw/moss.g` in BRL-CAD itself.
Set `FLTKCAD_TEST_DATABASE` to override it; database tests are skipped if absent.
Tcl scripting follows `BRLCAD_ENABLE_TCL`, and the legacy Tk lifecycle test
also requires `BRLCAD_ENABLE_TK`.

## Xmin testing

For deterministic headless testing, build FLTK and BRL-CAD against the same
Xmin SDK and configure with `BRLCAD_X11_PROVIDER=XMIN`.
CTest discovers `xmin-run` and `xminctl` from that SDK; their paths may be
overridden with `FLTKCAD_XMIN_RUN` and `FLTKCAD_XMINCTL`. The test suite exercises
both display backends, real Tk lifecycle, the terminal, HTML loading, editor
transactions, failure limits, and the reusable package boundary.

Detailed component behavior is documented in [docs/html.md](docs/html.md),
[docs/terminal.md](docs/terminal.md), and [docs/display.md](docs/display.md).
