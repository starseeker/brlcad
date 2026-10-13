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

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

#include <string>

extern "C" {
#include "dm.h"
}

int
main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  const auto backend = argc > 1 && std::string(argv[1]) == "--software"
                         ? Fl_Cad_View::Backend::OSMesa
                         : Fl_Cad_View::Backend::Auto;
  Fl_Double_Window window(800, 600, "FLTK libdm view");
  Fl_Cad_View view(0, 0, window.w(), window.h(), backend);
  view.draw_handler([](struct bview* bv) {
    auto* dmp = static_cast<struct dm*>(bv->dmp);
    dm_set_fg(dmp, 255, 255, 0, 1, 1.0);
    dm_draw_line_2d(dmp, -0.8, -0.8, 0.8, 0.8);
    dm_draw_line_2d(dmp, -0.8, 0.8, 0.8, -0.8);
  });
  window.resizable(view);
  window.end();
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
