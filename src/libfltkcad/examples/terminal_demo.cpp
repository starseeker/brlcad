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

#include <algorithm>
#include <string>
#include <vector>

int
main(int, char** argv)
{
  bu_setprogname(argv[0]);
  Fl_Double_Window window(800, 500, "FLTK CAD terminal");
  Fl_Cad_Terminal terminal(0, 0, window.w(), window.h());
  window.resizable(terminal);
  window.end();

  terminal.append("FLTK CAD terminal demo\n");
  terminal.command_handler([&terminal](const std::string& command) {
    terminal.append("command: " + command + "\n");
  });
  terminal.completion_handler([](const std::string& input, std::size_t cursor) {
    const std::vector<std::string> commands = {"clear", "close", "draw", "help"};
    const auto separator = input.rfind(' ', cursor == 0 ? 0 : cursor - 1);
    const std::size_t begin = separator == std::string::npos ? 0 : separator + 1;
    const std::string seed = input.substr(begin, cursor - begin);

    Fl_Cad_Terminal::Completion result{begin, cursor, {}};
    std::copy_if(commands.begin(), commands.end(),
                 std::back_inserter(result.candidates),
                 [&seed](const std::string& command) {
                   return command.compare(0, seed.size(), seed) == 0;
                 });
    return result;
  });

  window.show();
  return Fl::run();
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
