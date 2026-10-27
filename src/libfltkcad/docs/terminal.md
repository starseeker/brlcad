# Command terminal

`Fl_Cad_Terminal` is a command console built with standard `Fl_Text_Display`
and `Fl_Text_Buffer` widgets. The transcript, prompt, and pending command share
one text view. The transcript is protected from editing. Asynchronous output
preserves the pending command, insertion point, and selections in history.
Mouse selection, copying, and scrolling use FLTK's text widget behavior.

Enter submits; Up/Down and Control-P/N navigate history; Tab completes.
Left/Right, Home/End, Control-A/E/B/F, Backspace, and Delete edit the command.
Control-U/K/W erase the prefix, suffix, or preceding word. Control-L clears the
screen and Control-C cancels input and calls the optional interrupt handler.
Control-V, Shift-Insert, and the middle mouse button paste. Select output with
the mouse and use Control-Shift-C to copy. Shift navigation selects pending
input; pasting or deleting replaces that selection, and Control-X cuts it.
Pasted newlines stay in the pending
command until Enter submits it.

The widget deliberately does not depend on libged. Applications provide a
command callback and an optional completion callback:

```cpp
Fl_Cad_Terminal terminal(0, 0, 800, 240);
terminal.command_handler([](const std::string& command) {
  execute_command(command);
});
terminal.completion_handler([](const std::string& input, std::size_t cursor) {
  return complete_command(input, cursor);
});
terminal.append("ready\n");
```

Completion results identify the byte range to replace and supply candidate
strings. The widget applies a single candidate directly; multiple candidates
extend the common prefix and open a selectable list. History navigation keeps
and restores the user's current draft. History ranges may be consolidated for
multi-step commands such as libged's `GED_MORE` flow.

## Width and colors

Historical output and pending input wrap to the current viewport by default,
including long words and object names. Resizing changes only their display:
explicit newlines, table rows, and copied text retain their original content.
Widening the view reverses visual wrapping without rerunning commands.
`wrap_history(false)` disables wrapping and enables horizontal scrolling.
Fltkged exposes this as the persistent **Wrap historical output** preference.
`columns()` reports the current viewport width in character cells.

The widget defaults to black text on white. `colors(foreground, background)`
updates both new and existing text using the default palette. ANSI SGR colors
(base, bright, indexed, and RGB), bold, italic, underline, and strikeout are
translated into native text styles. Carriage return, backspace, and erase to
end of line support progress output. Escape sequences may span appends.
Other terminal cursor controls and OSC sequences are ignored: this widget is
a command transcript, rather than a full-screen terminal emulator. Up to 62
distinct text styles are retained; additional styles use the default palette.

No terminal wrapping patch or private FLTK API is required. `output()` returns
the native text view for focus and selection; `text()` returns a snapshot of
the transcript and live prompt/input with their original newlines.

Fltkged passes `columns()` to `ged_set_terminal_columns()` before executing
commands. Libged's default width remains 80 for callers that do not supply a
viewport. A terminal width is presentation state, independent of geometry,
selection, and command parsing.

## Failure and performance behavior

The default limits retain 1,000 submitted commands and 10,000 output lines.
Retained output is also bounded to 16 MiB, including unterminated lines.
Input is limited to 64 KiB, one `append()` call to 4 MiB, prompts to 4 KiB, and
completion results to 200 candidates. Oversized values are truncated or
rejected with a diagnostic available from `last_error()`; output receives an
explicit truncation marker. Limits count logical lines and bytes, so resizing
does not discard history.

Command, completion, interrupt, and error callback exceptions are caught so they do not
unwind through the FLTK event loop. Applications should still split sustained
high-volume output into event-loop-sized batches—the per-call limit prevents a
single accidental allocation, while the transcript line and byte limits bound
retained history.

## Build and test

The terminal is part of `libfltkcad`. Its history tests do not need a display;
widget tests run under Xmin when `BRLCAD_X11_PROVIDER=XMIN` is selected.
