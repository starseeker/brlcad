# Display manager view

`Fl_Cad_View` is a single FLTK geometry-view widget with two interchangeable
libdm backends:

- `SystemOpenGL` uses an `Fl_Gl_Window` and BRL-CAD's toolkit-neutral `extgl`
  display manager;
- `OSMesa` uses BRL-CAD's `swrast` display manager and copies its off-screen
  RGBA buffer into a normal FLTK widget; and
- `Auto` tries system OpenGL and defers a clean fallback to OSMesa if context
  or display-manager initialization fails.

Both paths share one render lifecycle, view-input handling, invalidation,
capture, and error surface. The widget owns a default `bview`; `view(borrowed)`
may replace it, but the caller must keep a borrowed view alive until it is
replaced or the widget is destroyed. The `dm` and framebuffer accessors are
non-owning.

```cpp
Fl_Cad_View view(0, 0, 800, 600, Fl_Cad_View::Backend::Auto);
view.error_handler([](const std::string& message) {
  log_error(message);
});
view.invalidate();
```

The built-in bindings follow MGED: left/right clicks zoom out/in by two;
Shift drags translate, Control drags rotate, and Shift-Control drags scale
with any of the three mouse buttons. Adding Alt constrains motion to X, Y, or
Z according to the button. The wheel uses MGED's zoom direction and factors.
Common letter/number shortcuts select standard views. They may be disabled when an application supplies its own edit
tools. A custom drawing callback replaces the standard `dm_draw_objs` scene
pass; libdm faceplate overlays are still drawn afterwards.

`bview::gv_s` controls the existing libdm view-parameter, scale, axes, grid,
ADC, and rectangle overlays. The reusable widget leaves their normal libbv
defaults intact; fltkged enables view parameters by default. The status font
size defaults to zero (automatic): libdm measures the text with the backend,
chooses a size from 10 to 14, and wraps at spaces when a single line would need
smaller text. A positive `gv_view_params.font_size` fixes the size while retaining
wrapping. `view faceplate params font_size auto` (or `0`) restores automatic
sizing; `view faceplate params font_size 12` selects a fixed size. A viewport must
still have room for the resulting readable lines.

`faceplate_gui(description)` registers optional menu and slider controls from
`dm/faceplate.h`. Their providers supply labels, colors, values, and selection
callbacks without depending on Tcl/Tk or librt. The description and callback
data are borrowed and must remain alive until `faceplate_gui(nullptr)` or
widget destruction. Registration survives backend recreation. Middle-button
events on visible controls are dispatched before application tools; slider
drags retain their starting row and clamp to its horizontal bounds. Fltkged
keeps application actions and knob state outside the renderer.

## BRL-CAD integration

The widget uses BRL-CAD's native external-context display manager in
`src/libdm/extgl`. Display plugins load with isolated symbol binding where the
platform supports it, so Xmin's GL and the swrast plugin can each use their
own OSMesa implementation.

The host toolkit and BRL-CAD `libdmgl` must use the same OpenGL ABI. For an
Xmin FLTK build, configure BRL-CAD with `BRLCAD_X11_PROVIDER=XMIN` and the same
`Xmin_DIR`.

## Failure and performance behavior

Explicit backend selection fails with a diagnostic instead of silently
changing renderers. `Auto` alone permits fallback. Context creation, plugin
loading, resize, frame begin/end, capture, and callback failures are surfaced
through `last_error()` and the error callback; application callback exceptions
are contained.

Framebuffer dimensions default to at most 4096 by 4096 and 16 million pixels.
The limits are checked before allocating or resizing. BRL-CAD's swrast backend
allocates an exact-size framebuffer and enforces OSMesa's dimension ceiling. OSMesa frames are cached
until invalidated. OpenGL captures render a fresh frame and read its back
buffer before FLTK presents it.

Headless tests render and capture a nonblank scene through both backends under
Xmin, exercise callback containment, and verify dimension-limit failures.
