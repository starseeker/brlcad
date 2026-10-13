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

#include <FL/Fl.H>

#include <iostream>

extern "C" {
#include "bu/app.h"
}

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  Fl_Ged_Tcl interpreter;
  if (!interpreter.available()) {
    std::cerr << interpreter.last_error() << '\n';
    return 1;
  }
  if (!interpreter.enable_legacy_tk() || !interpreter.legacy_tk_active()) {
    std::cerr << interpreter.last_error() << '\n';
    return 1;
  }
  if (interpreter.evaluate(
        "wm withdraw .; toplevel .fltkged_test; update idletasks; "
        "destroy .fltkged_test; expr {[winset] eq \"mged\" && "
        "[llength $mged_players] == 1}") != 0 ||
      interpreter.result() != "1") {
    std::cerr << interpreter.last_error() << '\n';
    return 1;
  }
  Fl::check();
  return 0;
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
