/* BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */

#include "common.h"
#include "bu/app.h"

#include "fltkcad/Fl_Cad_Terminal.H"

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Terminal.H>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void
require(bool condition, const char* message)
{
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

std::string transcript(const Fl_Cad_Terminal& terminal)
{
  const char* buffer = terminal.output().text(true);
  const std::string text = buffer ? buffer : "";
  std::free(const_cast<char*>(buffer));
  return text;
}

void event(Fl_Cad_Terminal& terminal, int type, int key = 0,
           const std::string& text = {}, int state = 0)
{
  Fl::e_number = type;
  Fl::e_keysym = key;
  Fl::e_original_keysym = key;
  Fl::e_text = const_cast<char*>(text.c_str());
  Fl::e_length = static_cast<int>(text.size());
  Fl::e_state = state;
  require(terminal.output().handle(type) != 0, "terminal did not handle editor event");
  Fl::e_text = nullptr;
  Fl::e_length = 0;
}

} // namespace

int
main(int, char** argv)
{
  bu_setprogname(argv[0]);
  if (std::getenv("DISPLAY") == nullptr) {
    return 77;
  }
  Fl_Double_Window window(640, 360);
  Fl_Cad_Terminal terminal(0, 0, window.w(), window.h());
  window.end();

  std::string submitted;
  terminal.command_handler([&submitted](const std::string& command) {
    submitted = command;
  });
  terminal.input("draw all");
  terminal.submit();

  require(submitted == "draw all", "command callback did not receive input");
  require(terminal.input().empty(), "input was not cleared after submission");
  require(terminal.history_size() == 1, "submitted command was not recorded");
  require(terminal.history_at(0) == "draw all", "recorded command is wrong");
  require(transcript(terminal).find("mged> draw all") != std::string::npos,
          "submitted command was not echoed");

  terminal.input("he");
  terminal.completion_handler([](const std::string&, std::size_t cursor) {
    return Fl_Cad_Terminal::Completion{0, cursor, {"help"}};
  });
  terminal.complete();
  require(terminal.input() == "help", "single completion was not applied");

  terminal.append("asynchronous output\n");
  require(terminal.input() == "help", "output disturbed pending input");

  for (int child = 0; child < terminal.children(); ++child) {
    require(dynamic_cast<Fl_Input*>(terminal.child(child)) == nullptr,
            "terminal still uses a separate input widget");
  }
  terminal.clear();
  terminal.output().clear_history();
  terminal.input("drw");
  event(terminal, FL_KEYDOWN, FL_Home);
  event(terminal, FL_KEYDOWN, FL_Right);
  event(terminal, FL_KEYDOWN, FL_Right);
  event(terminal, FL_PASTE, 0, "a");
  require(terminal.input() == "draw", "insertion did not use the terminal cursor");
  require(transcript(terminal).find("mged> draw") != std::string::npos,
          "live command is missing from the terminal display");
  std::size_t position = 0;
  terminal.completion_handler([&position](const std::string&, std::size_t cursor) {
    position = cursor;
    return Fl_Cad_Terminal::Completion{};
  });
  terminal.append("one ");
  terminal.append("two\n");
  terminal.complete();
  require(position == 3, "asynchronous output moved the insertion point");
  const std::string output = transcript(terminal);
  require(output.find("one two") != std::string::npos,
          "partial output was split by the live prompt");
  require(output.find("mged> draw") == output.rfind("mged> draw"),
          "editing duplicated the live command in scrollback");
  event(terminal, FL_KEYDOWN, FL_End);
  event(terminal, FL_PASTE, 0, " café");
  event(terminal, FL_KEYDOWN, FL_BackSpace);
  require(terminal.input() == "draw caf", "backspace split a UTF-8 character");
  event(terminal, FL_KEYDOWN, 'w', {}, FL_CTRL);
  require(terminal.input() == "draw ", "Control-W did not erase the preceding word");
  event(terminal, FL_KEYDOWN, 'u', {}, FL_CTRL);
  require(terminal.input().empty(), "Control-U did not clear the command prefix");
  event(terminal, FL_PASTE, 0, "hello");
  event(terminal, FL_KEYDOWN, 'a', {}, FL_CTRL);
  event(terminal, FL_KEYDOWN, FL_Delete);
  require(terminal.input() == "ello", "Delete did not erase at the cursor");
  event(terminal, FL_KEYDOWN, FL_Enter);
  require(submitted == "ello", "Enter did not submit from the terminal");
  event(terminal, FL_KEYDOWN, FL_Up);
  require(terminal.input() == "ello", "Up did not recall history in the terminal");
  event(terminal, FL_KEYDOWN, FL_Down);
  require(terminal.input().empty(), "Down did not restore the pending command");
  terminal.input(std::string(2000, 'x'));
  event(terminal, FL_KEYDOWN, FL_Home);
  event(terminal, FL_PASTE, 0, "a");
  terminal.resize(0, 0, 360, 180);
  event(terminal, FL_KEYDOWN, FL_End);
  terminal.append("background output\n");
  require(terminal.input() == "a" + std::string(2000, 'x'),
          "wrapping, resize, or background output damaged a long command");
  event(terminal, FL_KEYDOWN, 'c', {}, FL_CTRL);
  require(terminal.input().empty(), "Control-C did not cancel pending input");

  terminal.clear();
  terminal.output().clear_history();
  terminal.resize(0, 0, 800, 180);
  const int wide_columns = terminal.columns();
  const std::string wide_line(static_cast<std::size_t>(wide_columns - 1), 'W');
  terminal.append(wide_line + "\n");
  terminal.input(std::string(200, 'i'));
  terminal.resize(0, 0, 240, 180);
  const int narrow_columns = terminal.columns();
  require(narrow_columns < wide_columns, "shrink did not change active columns");
  require(terminal.output().cursor_col() < narrow_columns,
          "input cursor is outside the resized viewport");
  require(transcript(terminal).find(wide_line) != std::string::npos,
          "shrinking truncated wide history");
  terminal.input("");
  terminal.append(std::string(static_cast<std::size_t>(narrow_columns + 3), 'n'));
  const std::string narrowed = transcript(terminal);
  require(narrowed.find(std::string(static_cast<std::size_t>(narrow_columns), 'n') +
                         "\nnnn") != std::string::npos,
          "new output did not wrap at the visible column width");
  terminal.resize(0, 0, 960, 180);
  require(terminal.columns() > wide_columns &&
            terminal.output().display_columns() >= terminal.columns(),
          "width-only enlargement did not allocate usable columns");
  terminal.input(std::string(200, 'i'));
  terminal.resize(0, 0, 36, 180);
  require(terminal.columns() >= 1 &&
            terminal.output().cursor_col() < terminal.columns(),
          "very narrow terminal left the editor cursor offscreen");
  terminal.resize(0, 0, 640, 360);
  terminal.input("");
  require(terminal.output().textfgcolor_default() == FL_BLACK &&
            terminal.output().textbgcolor_default() == FL_WHITE,
          "terminal does not default to black on white");
  terminal.output().append("\x1b[31m");
  const Fl_Color ansi_foreground = terminal.output().textfgcolor();
  terminal.colors(FL_WHITE, FL_BLACK);
  require(terminal.output().textfgcolor() == ansi_foreground,
          "changing the palette interrupted an active ANSI color");
  terminal.output().append("\x1b[0m");
  require(terminal.output().textfgcolor_default() == FL_WHITE &&
            terminal.output().textbgcolor_default() == FL_BLACK &&
            transcript(terminal).find(wide_line) != std::string::npos,
          "changing terminal colors damaged output or failed to change defaults");
  terminal.colors(FL_BLACK, FL_WHITE);

  terminal.command_handler([](const std::string&) {
    throw std::runtime_error("intentional command failure");
  });
  terminal.input("throw");
  terminal.submit();
  require(terminal.last_error().find("command callback failed") !=
            std::string::npos,
          "command callback exception escaped or was not reported");
  terminal.completion_handler([](const std::string&, std::size_t) {
    throw std::runtime_error("intentional completion failure");
    return Fl_Cad_Terminal::Completion{};
  });
  terminal.input("x");
  terminal.complete();
  require(terminal.last_error().find("completion callback failed") !=
            std::string::npos,
          "completion callback exception escaped or was not reported");

  auto limits = terminal.limits();
  limits.max_input_bytes = 4;
  limits.max_append_bytes = 4;
  terminal.limits(limits);
  terminal.input("oversized");
  require(terminal.input() == "over", "oversized input was not bounded");
  require(!terminal.last_error().empty(), "input truncation was not reported");
  terminal.append("oversized");
  require(terminal.last_error().find("append") != std::string::npos,
          "output truncation was not reported");

  limits.max_prompt_bytes = 4;
  terminal.limits(limits);
  terminal.error_handler([](const std::string&) {
    throw std::runtime_error("intentional error-handler failure");
  });
  terminal.prompt("oversized prompt");
  require(terminal.prompt() == "over", "oversized prompt was not bounded");

  terminal.prompt("four");
  limits.max_prompt_bytes = 2;
  terminal.limits(limits);
  require(terminal.prompt() == "fo",
          "lowering the prompt limit did not bound the current prompt");

  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
