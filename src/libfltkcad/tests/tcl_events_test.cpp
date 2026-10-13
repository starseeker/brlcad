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
#include "fltkcad/Fl_Ged_Tcl.H"
#include "bu/app.h"
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <chrono>
#include <iostream>
#include <string>

namespace {
bool require(bool condition, const std::string& message)
{
  if (!condition) std::cerr << message << '\n';
  return condition;
}
}

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  Fl_Double_Window window(200, 100);
  window.show();
  Fl_Ged_Tcl interpreter;
  if (!require(interpreter.available(), interpreter.last_error())) return 1;
  interpreter.event_processing(true);
  std::string output;
  interpreter.output_handler([&](const std::string& text, bool) { output += text; });
  if (!require(interpreter.evaluate("puts -nonewline {hello }; puts world; puts stderr diagnostic") == 0 &&
               output == "hello world\ndiagnostic\n", "Tcl puts did not reach the output handler")) return 1;
  if (!require(!interpreter.command_complete("set value {first") &&
               interpreter.command_complete("set value {first\nsecond}"),
               "Tcl multiline completeness is incorrect")) return 1;

  bool fltk_timer = false;
  Fl::add_timeout(0.01, [](void* data) { *static_cast<bool*>(data) = true; }, &fltk_timer);
  if (!require(interpreter.evaluate("after 50 {set done 1}; vwait done; update") == 0 && fltk_timer,
               "FLTK timers did not run while Tcl waited: " + interpreter.last_error())) return 1;
  bool tcl_timer = false;
  interpreter.ui_handler([&](const std::vector<std::string>&, std::string&) {
    tcl_timer = true;
    return 0;
  });
  if (interpreter.evaluate("after 10 {fltkged tick}") != 0) return 1;
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
  while (!tcl_timer && std::chrono::steady_clock::now() < deadline) Fl::wait(0.01);
  if (!require(tcl_timer, "Tcl timers did not run in FLTK's event loop")) return 1;
  tcl_timer = false;
  fltk_timer = false;
  Fl::add_timeout(0.01, [](void* data) { *static_cast<bool*>(data) = true; }, &fltk_timer);
  if (interpreter.evaluate("after 0 {after 50 {set nested 1}; vwait nested; fltkged tick}") != 0) return 1;
  const auto nested_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
  while (!tcl_timer && std::chrono::steady_clock::now() < nested_deadline) Fl::wait(0.01);
  if (!require(tcl_timer && fltk_timer, "A Tcl timer's vwait starved FLTK events")) return 1;
  Fl::add_timeout(0.01, [](void* data) { static_cast<Fl_Ged_Tcl*>(data)->interrupt(); }, &interpreter);
  if (!require(interpreter.evaluate("vwait never_set") != 0 &&
               interpreter.evaluate("expr {1 + 1}") == 0 && interpreter.result() == "2",
               "Interrupting vwait prevented later commands from running")) return 1;
  if (!require(interpreter.evaluate(R"tcl(
    lassign [chan pipe] reader writer
    fconfigure $reader -blocking 0
    fileevent $reader readable {set message [gets $reader]; fileevent $reader readable {}}
    after 10 {puts $writer ready; flush $writer}
    vwait message
    close $reader
    close $writer
    expr {$message eq "ready"}
  )tcl") == 0 && interpreter.result() == "1",
    "Tcl file events did not run: " + interpreter.last_error())) return 1;
  interpreter.event_processing(false);
  for (int script = 0; script < 3; ++script) {
    if (!require(interpreter.evaluate("for {set i 0} {$i < 200000} {incr i} {set value $i}") == 0,
                 "Tcl command budget accumulated across evaluations: " + interpreter.last_error())) return 1;
  }
  if (!require(!interpreter.legacy_tk_active(), "Tcl event handling unexpectedly loaded Tk")) return 1;
  window.hide();
  return 0;
}
