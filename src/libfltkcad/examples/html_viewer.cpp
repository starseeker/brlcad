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

#include <fltkcad/Fl_Litehtml_View.H>

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

#include <iostream>

int main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  Fl_Double_Window window(900, 700, "FLTK litehtml viewer");
  Fl_Litehtml_View view(0, 0, window.w(), window.h());
  window.resizable(view);
  window.end();

  view.error_handler([](const std::string& message) {
    std::cerr << message << '\n';
  });
  view.link_handler([](const std::string& uri) {
    std::cout << "Activated link: " << uri << '\n';
  });

  if(argc > 1)
  {
    if(!view.load(argv[1]))
    {
      return 1;
    }
    if(argc > 2 && view.find(argv[2]) == 0)
    {
      std::cerr << "No matches for: " << argv[2] << '\n';
      return 2;
    }
  }
  else
  {
    view.html("<h1>FLTK + litehtml</h1>"
              "<p>A native FLTK HTML widget with <a href='#search'>links</a>, "
              "images, scrolling, and document search.</p>"
              "<h2 id='search'>Search</h2><p>Try the public find API.</p>");
  }

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
