# TkLiteHtml

TkLiteHtml is a native Tcl/Tk widget backed by
[litehtml](https://github.com/litehtml/litehtml). It is intended to
provide the HTML viewing, link handling, search, find, and highlighting
capabilities needed to replace tkhtml3 and Hv3 in BRL-CAD.

The renderer deliberately has a small dependency surface: Tcl/Tk and
litehtml. It does not use Cairo, Pango, or a platform graphics API
directly. Backgrounds, images, borders, gradients, and highlights are
composed into an RGBA buffer and uploaded through Tk's photo API; text
measurement and drawing use Tk's public font APIs. A canvas image paints
the background and text together, letting Tk buffer each complete frame
through its portable API. This keeps the implementation applicable to
X11, Windows, macOS Aqua, and other Tk ports.

## Build

Within BRL-CAD, the native widget and Tcl browser helpers are built into
libtclcad when Tcl/Tk is enabled. `tclcad_init(interp, 1, ...)` initializes
TkLiteHtml in that interpreter. No separate widget module or Tcl script
is needed at runtime. Headless initialization does not initialize Tk.

The widget supports litehtml v0.10, supplied by bext, and the newer
API with redraw callbacks. CMake detects which API is available when
configuring. BRL-CAD discovers litehtml through its usual bext package
handling and `FindLITEHTML.cmake`, which uses the installed CMake package
and its Gumbo dependency when available. Installations without package
configuration files are also supported.

This directory can also be copied and built independently. It only
requires Tcl/Tk 8.6 and an installed litehtml CMake package:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/dependencies
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /desired/prefix
```

The standalone build installs a loadable Tcl package. Both builds embed
`library/tklitehtml.tcl` into the native package, so `package require
tklitehtml 0.1.0` provides the same widget and browser API. Changes to the
Tcl script cause CMake to regenerate the embedded copy.

The standalone build uses Tcl and Tk stubs. If CMake cannot find the intended
Tcl/Tk installation, set `TCL_LIBRARY`, `TK_LIBRARY`,
`TCL_STUB_LIBRARY`, `TK_STUB_LIBRARY`, `TCL_INCLUDE_PATH`, and
`TK_INCLUDE_PATH` explicitly. Set `TK_PLATFORM_INCLUDE_PATH` when
platform headers referenced by `tk.h` are outside the compiler's normal
include path. These are build inputs only; TkLiteHtml does not link
directly to X11 or Xmin.

GUI tests require a display. They return CTest's skip code when Tk
cannot open one. The same suite can be run under Xvfb, Xmin, or a
native desktop session.

## Basic use

```tcl
package require Tk
package require tklitehtml 0.1.0

scrollbar .scroll -command {.html yview}
tklitehtml::widget .html \
    -yscrollcommand {.scroll set} \
    -linkcommand openLink \
    -resourcecommand loadResource

pack .scroll -side right -fill y
pack .html -side left -fill both -expand 1

.html load -baseurl https://example.invalid/help/index.html {
    <h1>Help</h1>
    <p>Read the <a href="commands.html">command reference</a>.</p>
}
```

Mouse-wheel bindings are installed for the `TkLiteHtml` widget class.
They normalize X11 button events and Win32/Aqua `<MouseWheel>` deltas,
including Shift-wheel horizontal scrolling, into the widget's `xview`
and `yview` commands.

## Widget API

Create a widget with:

```tcl
tklitehtml::widget pathName ?option value ...?
```

Supported options are:

- `-width`, `-height`, and `-background`
- `-maxdocumentbytes`, `-maximagepixels`, and `-maxviewportpixels`
- `-xscrollcommand` and `-yscrollcommand`
- `-linkcommand`
- `-resourcecommand`

The three resource limits default to 16 MiB, 16 million pixels, and
16 million pixels, respectively. They bound the combined document and
local resource input, the combined decoded images, and the backing
image used for the visible viewport. A failed `load` leaves the prior
document intact.

The widget command supports:

```text
pathName load ?-baseurl uri? ?-userstyle css? html
pathName clear
pathName configure ?option? ?value ...?
pathName cget option
pathName text
pathName search ?-all? ?-nocase? ?--? pattern
pathName highlight configure name -background color
pathName highlight set name ranges
pathName highlight clear ?name?
pathName find set ?-case? pattern
pathName find current|next|previous|clear
pathName fragment identifier
pathName bbox index
pathName see index
pathName index index
pathName selection set first last
pathName selection get|range|clear
pathName title
pathName xview ?moveto fraction|scroll count units|pages?
pathName yview ?moveto fraction|scroll count units|pages?
```

Text indexes are Unicode character offsets. Search and highlight ranges
are half-open `{first last}` pairs, so `last` is the first character
after the match. Text follows document order, with structural line breaks
and CSS whitespace handling. Soft wrapping adds no characters, so text
and ranges remain stable when the viewport is resized. Character boxes
are rebuilt for the current layout. `search` returns the first range unless
`-all` is specified. `find set` highlights all matches, selects the first, and
returns the match count; `find next` and `find previous` cycle the
current match and scroll it into view. `fragment` scrolls to an element
by HTML `id`, and `title` returns the current document title.

An index may be a character offset, `end`, or `@x,y` in viewport pixels.
Coordinate indexes resolve to the nearest text insertion position and
account for scrolling. `selection set` accepts either direction;
`selection range` returns a half-open `{first last}` pair or an empty
list, and `selection get` returns the selected text. Selection is
independent of find matches and survives reflow. Loading another document
or calling `clear` clears selection; a failed load preserves it.

Class bindings provide drag selection, double-click word selection,
triple-click paragraph selection, Shift-click extension, and scrolling
when dragging outside the viewport. `<<Copy>>` copies to CLIPBOARD, and
Control-A (Command-A on Aqua) selects all. On X11, selecting also exports
PRIMARY for middle-button pasting in other applications. Selecting alone
does not replace CLIPBOARD. Dragging over links suppresses navigation.

`<<Selection>>` reports selection changes; `<<DocumentChanged>>` reports
successful document replacement or clearing. The `selection`, `find`,
and `find-current` highlight names are used by the built-in operations.

## Local document browser

For applications that need navigation rather than a single rendered
document, the package also supplies a composite browser:

```tcl
set browser [tklitehtml::browser .browser]
pack $browser -fill both -expand 1
$browser goto /path/to/index.html
```

The browser adds scrollbars and supports `goto`, `back`, `location`,
`title`, `error`, `findbar ?pattern?`, `autofind boolean ?pattern?`, and
`html` (which returns the underlying widget).
Relative and fragment links, file URIs, percent-encoded paths, and
navigation history are handled without an HTTP client. A
`-linkcommand` callback may replace the default navigation; it receives
the browser path, resolved URL, and link attribute dictionary. The
portable `tklitehtml::file_uri` and `tklitehtml::local_path` commands
convert between local paths and file URIs.

Right-click the HTML view for Copy, Select All, Find, Find Next, Find
Previous, and Clear Highlights. On Aqua, Control-click also opens the
menu. Copy is enabled when text is selected.

Control-F (Command-F on Aqua) opens the find bar, using selected text as
the initial query when available. `findbar ?pattern?` opens the same bar
programmatically and returns its entry widget. Searches update while
typing, highlight all matches, and show the current match count. Match
case controls case sensitivity. Next/Previous, F3/Shift-F3, and
Return/Shift-Return in the entry cycle through matches. Escape or Close
hides the bar and clears find highlights. Clear Highlights clears the
query without changing text selection. An open query is reapplied when
another document is loaded successfully; failed navigation preserves
the current results.

`autofind true ?pattern?` lets an embedding application show and update
the find bar without moving keyboard focus. Closing an automatically opened
bar suppresses it until `autofind false` resets the automatic-find cycle.
An explicit `findbar` request remains open when automatic find is disabled.

Navigation reads local files in bounded chunks. Missing, unreadable,
oversized, or invalid pages leave the current page and Back history
intact. The browser reports the error in a visible banner and returns a
Tcl error; successful navigation clears the banner. Calling `error`
returns the current message, and `error message` lets an embedding
application report its own load failure in the same banner.

## Local help viewer

Applications can build a standard documentation view around the browser:

```tcl
set help [tklitehtml::helpviewer .help \
    -home /path/to/start.html \
    -contents /path/to/contents.html]
pack $help -fill both -expand true
```

The toolbar provides Back, Home, Find, and Contents controls. When a
contents document is supplied, it appears beside the main document in a
paned window. The sash adjusts their relative widths, and Contents toggles
the sidebar without discarding its document or width. Links selected in the
contents pane open in the main pane.

The widget command supports `back`, `home`, `goto`, `findbar`, and
`showcontents ?boolean?`. `browser` and `contentsbrowser` return the
underlying browser widgets. Other browser operations are forwarded to the
main document browser. Only local files and fragment links are supported;
network retrieval is intentionally outside this component's scope.

### Link callback

The command prefix supplied with `-linkcommand` receives the widget
path, resolved URL, and a dictionary of available `href`, `title`, `id`,
and `class` attributes:

```tcl
proc openLink {widget url attributes} {
    puts "open $url"
}
```

Fragment links scroll within the current document before invoking the
callback.

### Resource callback

The command prefix supplied with `-resourcecommand` receives the widget
path, resource kind, and resolved URL:

```tcl
proc loadResource {widget kind url} {
    switch -- $kind {
        image      { return myTkPhotoImage }
        stylesheet { return {body { color: #222; }} }
    }
    return {}
}
```

For `image`, return the name of a Tk photo image. For `stylesheet`,
return CSS text. When the callback returns nothing, local image and
stylesheet paths are loaded directly. Network access is intentionally
left to the embedding application.

## Current rendering scope

The widget implements the core static-document path: Tk fonts and text,
software-composited colors and transparency, rectangular borders, list
markers, Tk photo images, linear gradients, clipping, links, scrolling,
and responsive relayout. Radial and conic gradients currently use a
solid-color fallback, and rounded border rasterization is not yet
implemented. JavaScript, forms, and an HTTP client are outside the
widget's scope.
