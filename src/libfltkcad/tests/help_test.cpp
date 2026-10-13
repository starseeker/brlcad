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
#include "html_test_utils.h"
#include "fltkcad/Fl_Ged_Help.H"
#include "bu/app.h"
#include "bu/file.h"
#include <FL/Fl.H>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
void require(bool condition, const char* message)
{
  if (!condition) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
}

int main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  if (!std::getenv("DISPLAY")) return 77;
  require(argc == 3, "two HTML fixtures are required");
  Fl_Ged_Help viewer;
  require(viewer.load(argv[1]), "first HTML page did not load");
  const std::string first_title = viewer.label();
  require(viewer.load(argv[2]), "second HTML page did not load");
  const std::string second_title = viewer.label();
  require(first_title != second_title, "HTML fixture titles must differ");
  require(viewer.back() && first_title == viewer.label(), "Back did not restore the first page");
  require(!viewer.back(), "Back moved past the beginning of history");
  require(viewer.forward() && second_title == viewer.label(), "Forward did not restore the second page");
  require(!viewer.load("missing-manual.html") && !viewer.last_error().empty(),
          "missing manual was not diagnosed");
  require(viewer.back() && first_title == viewer.label(), "failed load damaged navigation history");
  require(viewer.load(argv[1]) && !viewer.forward(), "new navigation retained stale forward history");
  Fl_Litehtml_View* document = nullptr;
  for (int index = 0; index < viewer.children(); ++index) {
    document = dynamic_cast<Fl_Litehtml_View*>(viewer.child(index));
    if (document) break;
  }
  require(document, "manual HTML widget missing");
  const std::string base_uri = document->base_uri();
  require(document->html(
    "<title>Fragment navigation</title><style>body{margin:0}</style>"
    "<a href='#target'>needle</a><div style='height:700px'></div>"
    "<h2 id='target'>Another needle</h2>", base_uri), "navigation fixture did not load");
  require(document->find("needle") == 2, "navigation fixture search failed");
  click_first_html_link(*document);
  require(document->document_title() == "Fragment navigation" && document->match_count() == 2,
          "same-document hyperlink reloaded the page or lost its search state");
  require(!viewer.forward(), "fragment navigation changed page history");
  const std::string manual = bu_dir(nullptr, 0, BU_DIR_DOC, "html", "man1", "mged.html", nullptr);
  if (bu_file_exists(manual.c_str(), nullptr)) {
    require(viewer.manuals("en") && viewer.manual_count() > 0,
            "installed manual pages were not indexed");
    require(viewer.manuals("C.", 'n'), "C locale did not fall back to English manuals");
    require(!viewer.manuals("en", '2'), "unsupported manual section was accepted");
    require(viewer.load_default(), "default MGED documentation did not load");
  }
  viewer.show();
  Fl::check();
  viewer.resize(viewer.x(), viewer.y(), 1100, 800);
  Fl::check();
  viewer.hide();
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
