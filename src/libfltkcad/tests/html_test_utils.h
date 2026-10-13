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

#ifndef FLTKCAD_HTML_TEST_UTILS_H
#define FLTKCAD_HTML_TEST_UTILS_H

#include "fltkcad/Fl_Litehtml_View.H"
#include <FL/Fl.H>
#include <initializer_list>

// Navigation fixtures place their first link at the document origin.
inline void click_first_html_link(Fl_Litehtml_View& view)
{
  constexpr int link_inset = 10;
  Fl::e_x = view.x() + link_inset;
  Fl::e_y = view.y() + link_inset;
  Fl::e_keysym = FL_Button + FL_LEFT_MOUSE;
  for (int event : {FL_MOVE, FL_PUSH, FL_RELEASE}) {
    Fl::e_number = event;
    static_cast<Fl_Widget&>(view).handle(event);
  }
}

#endif

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
