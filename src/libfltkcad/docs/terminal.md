# Command terminal

`Fl_Cad_Terminal` edits commands directly in `Fl_Terminal`. The prompt, pending
command, cursor, and output share one terminal surface. Asynchronous output
restores the pending command and insertion point, including wrapped input.
Very long commands display the portion containing the cursor. Submitted
commands and output use FLTK's bounded scrollback.

Enter submits; Up/Down and Control-P/N navigate history; Tab completes.
Left/Right, Home/End, Control-A/E/B/F, Backspace, and Delete edit the command.
Control-U/K/W erase the prefix, suffix, or preceding word. Control-L clears the
screen and Control-C cancels input and calls the optional interrupt handler.
Control-V, Shift-Insert, and the middle mouse button paste. Select output with
the mouse and use Control-Shift-C to copy. Pasted newlines stay in the pending
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

Input and new output wrap at the current viewport width on every resize,
including width changes that leave the height unchanged. `columns()` reports
that width in character cells. FLTK's ring buffer retains its widest allocation
so resizing never truncates older output. Old rows retain their original line
breaks and can be read with the horizontal scrollbar; resizing does not rerun
commands or reformat tables already printed.

The widget defaults to black text on white. `colors(foreground, background)`
changes the default palette and updates existing text using the old defaults.
ANSI color and cursor sequences continue to be handled by `Fl_Terminal`.
The bext FLTK wrapping patch supplies an output margin independent of history
storage; an unpatched FLTK package does not provide this API.

Fltkged passes `columns()` to `ged_set_terminal_columns()` before executing
commands. Libged's default width remains 80 for callers that do not supply a
viewport. A terminal width is presentation state, independent of geometry,
selection, and command parsing.

## Failure and performance behavior

The default limits retain 1,000 submitted commands and 10,000 output lines.
Input is limited to 64 KiB, one `append()` call to 4 MiB, prompts to 4 KiB, and
completion results to 200 candidates. Oversized values are truncated or
rejected with a diagnostic available from `last_error()`; output receives an
explicit truncation marker. Appends are passed to FLTK in bounded chunks.

Command, completion, interrupt, and error callback exceptions are caught so they do not
unwind through the FLTK event loop. Applications should still split sustained
high-volume output into event-loop-sized batches—the per-call limit prevents a
single accidental allocation, while `Fl_Terminal`'s scrollback limit bounds
retained history.

## Build and test

The terminal is part of `libfltkcad`. Its history tests do not need a display;
widget tests run under Xmin when `BRLCAD_X11_PROVIDER=XMIN` is selected.
