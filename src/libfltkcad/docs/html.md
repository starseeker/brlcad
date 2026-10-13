# HTML view

`Fl_Litehtml_View` is a native FLTK HTML-view widget backed by
[litehtml](https://github.com/litehtml/litehtml). It is intended for portable,
embedded documentation viewers, including a possible replacement for the
tkhtml3/hv3 path used by BRL-CAD.

The public widget is `Fl_Litehtml_View`. The current implementation provides:

- HTML from a string or local file, with relative CSS and image resources;
- native FLTK text, backgrounds, gradients, borders, list markers, and images;
- vertical and horizontal scrolling with mouse and keyboard input;
- hover and link activation, with application-owned link policy;
- fragment navigation;
- rendered-text search with all matches highlighted and a distinct current
  match;
- callbacks for custom resource loading and error reporting; and
- headless interaction tests using Xmin.

The widget deliberately does not fetch network resources, launch external
applications, execute scripts, or implement application navigation. A host may
provide those policies through callbacks.

## Build

The HTML viewer is part of BRL-CAD's `libfltkcad`. From the BRL-CAD source root:

```sh
cmake -S . -B build \
  -DBRLCAD_EXT_DIR=/path/to/bext/build \
  -DBRLCAD_ENABLE_FLTK=ON \
  -DFLTKCAD_BUILD_EXAMPLES=ON
cmake --build build
```

FLTK and litehtml are discovered through their installed CMake packages in
bext. `BUILD_TESTING` controls the widget tests. Consumers use
`find_package(BRLCAD CONFIG REQUIRED COMPONENTS fltkcad)` and link
`BRLCAD::fltkcad`.

The viewer example accepts a file and an optional search string:

```sh
build/src/libfltkcad/examples/fltkcad_html_demo page.html optional-search
```

## Use

```cpp
#include <fltkcad/Fl_Litehtml_View.H>

Fl_Litehtml_View view(0, 0, 800, 600);
view.error_handler([](const std::string& message) {
  log_error(message);
});
view.link_handler([&view](const std::string& uri) {
  navigate_or_open(uri);
});

view.load("manual/index.html");
const std::size_t matches = view.find("rendering");
view.find_next();
```

`current_match()` is one-based and returns zero when no match is selected.
Search is byte-exact for UTF-8 text. Case-insensitive mode, the default, folds
ASCII letters without changing UTF-8 byte boundaries.

The default resource loader accepts local paths and `file://` URIs. A custom
loader can supply in-memory data for another URI scheme:

```cpp
view.resource_loader([](const std::string& uri)
  -> std::optional<Fl_Litehtml_View::Resource> {
  return load_from_application_resources(uri);
});
```

PNG, JPEG, GIF, SVG, and BMP images are decoded from loader results. CSS
imports use the same callback.

The document canvas defaults to white; applications may change the widget's
`color()`, and document CSS backgrounds are drawn above that canvas. Link
handlers run after litehtml finishes processing the click, so a handler can
load a new document safely.

## Failure and resource limits

`load()` accepts regular `.html`, `.htm`, and `.xhtml` files. A missing file,
unsupported extension, short read, or oversized file is rejected without
replacing the current document. `html()` rejects embedded NUL bytes and reports
parse or layout failures through `last_error()` and the optional error callback.
Exceptions from resource, link, and error callbacks are contained at the widget
boundary.

The `Limits` structure bounds document and resource bytes, resource requests,
cached and decoded images, extracted text, search work, fragment traversal, and
layout extents. Defaults include an 8 MiB document limit, 16 MiB per resource,
256 resources, 32 million decoded image pixels, and a one-million-pixel layout
extent. Limits can be adjusted before loading content:

```cpp
auto limits = view.limits();
limits.max_document_bytes = 2U * 1024U * 1024U;
limits.max_cached_images = 32;
view.limits(limits);
```

Resource loading is synchronous. Applications serving slow or remote content
should fetch and validate it outside the FLTK event thread, then expose cached
bytes through the resource callback.

## Headless tests with Xmin

Use an FLTK package built with Xmin support and the same Xmin SDK as BRL-CAD.
Configure from the BRL-CAD source root:

```sh
cmake -S . -B build-xmin \
  -DBRLCAD_EXT_DIR=/path/to/bext/build \
  -DBRLCAD_ENABLE_FLTK=ON \
  -DBRLCAD_X11_PROVIDER=XMIN \
  -DXmin_DIR=/path/to/xmin/lib/cmake/Xmin \
  -DFLTKCAD_BUILD_EXAMPLES=ON \
  -DBUILD_TESTING=ON
cmake --build build-xmin
ctest --test-dir build-xmin -R '^fltkcad_' --output-on-failure
```

When `xmin-run` and `xminctl` are alongside the Xmin server, CTest discovers
them automatically. The integration test launches the example in a real
headless FLTK session, loads external CSS and SVG resources, verifies the
current-search highlight in a captured frame, activates a fragment link,
scrolls through injected input, and confirms that each resulting frame changes.
It also audits the executable for accidental host X11/font-library linkage.
Normal builds do not require Xmin.

## Current scope

This is a focused documentation viewer rather than a browser engine. Rounded
background clipping and the more specialized CSS border styles are currently
approximated with FLTK primitives. The resource callback is synchronous, and
document-wide search belongs to the widget while cross-page indexing belongs
to the host application.
