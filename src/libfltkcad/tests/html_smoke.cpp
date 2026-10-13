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
#include "bu/app.h"

#include <fltkcad/Fl_Litehtml_View.H>

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{

void close_window(void* data)
{
  static_cast<Fl_Double_Window*>(data)->hide();
}

} // namespace

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  if(std::getenv("DISPLAY") == nullptr)
  {
    return 77;
  }

  Fl_Double_Window window(480, 320, "fltkcad-html-smoke");
  Fl_Litehtml_View view(0, 0, window.w(), window.h());
  window.resizable(view);
  window.end();
  if(view.color() != FL_WHITE) return 1;

  int resource_requests = 0;
  view.resource_loader([&resource_requests](const std::string& uri)
                         -> std::optional<Fl_Litehtml_View::Resource> {
    if(uri != "memory://manual/smoke.css")
    {
      return std::nullopt;
    }
    ++resource_requests;
    const std::string css = "body{font:16px sans-serif} mark{background:#ff0}";
    Fl_Litehtml_View::Resource resource;
    resource.data.assign(css.begin(), css.end());
    resource.mime_type = "text/css";
    return resource;
  });

  const std::string html =
    "<title>Smoke title</title><link rel='stylesheet' href='smoke.css'>"
    "<h1 id='top'>Smoke test</h1><p>Needle and another needle</p>"
    "<div style='height:700px'></div><h2 id='target'>Target</h2>";
  if(!view.html(html, "memory://manual/index.html") || resource_requests != 1 ||
     view.document_title() != "Smoke title" ||
     view.plain_text().find("Smoke test") == std::string::npos ||
     view.find("needle") != 2 || view.current_match() != 1 ||
     !view.find_next() || view.current_match() != 2 ||
     !view.find_previous() || view.current_match() != 1 ||
     !view.scroll_to_fragment("target"))
  {
    return 1;
  }
  view.clear_find();
  if(view.match_count() != 0 || view.current_match() != 0)
  {
    return 1;
  }

  // A browser's link handler replaces the document during navigation.
  const std::string linked_document =
    "<title>Next page</title><style>body{margin:0}</style>"
    "<a href='#target'>OPTIONS</a><div style='height:700px'></div>"
    "<h2 id='target'>Options</h2>";
  int clicked = 0;
  view.link_handler([&](const std::string&) {
    ++clicked;
    if(!view.html(linked_document, "memory://manual/next.html"))
      throw std::runtime_error("linked document did not load");
  });
  if(!view.html("<style>body{margin:0}</style><a href='next.html'>Next</a>",
                 "memory://manual/index.html")) return 1;
  click_first_html_link(view);
  if(clicked != 1 || view.document_title() != "Next page") {
    std::cerr << "Document replacement from a hyperlink failed\n";
    return 1;
  }
  click_first_html_link(view);
  if(clicked != 2 || view.document_title() != "Next page") {
    std::cerr << "Fragment hyperlink navigation failed\n";
    return 1;
  }
  view.link_handler({});

  if(view.load("not-html.txt") ||
     view.last_error().find("Unsupported HTML file type") == std::string::npos ||
     view.document_title() != "Next page")
  {
    return 1;
  }
  if(view.load("missing-document.html") ||
     view.last_error().find("HTML file") == std::string::npos ||
     view.document_title() != "Next page")
  {
    return 1;
  }

  view.resource_loader([](const std::string&)
                         -> std::optional<Fl_Litehtml_View::Resource> {
    throw std::runtime_error("intentional loader failure");
  });
  if(!view.html("<link rel='stylesheet' href='failure.css'><p>safe</p>",
                "memory://failure/index.html") ||
     view.last_error().find("Resource loader failed") == std::string::npos)
  {
    return 1;
  }

  auto search_limits = view.limits();
  search_limits.max_find_query_bytes = 2;
  view.limits(search_limits);
  if(view.find("oversized query") != 0 ||
     view.last_error().find("Find query exceeds") == std::string::npos)
  {
    return 1;
  }

  auto limits = view.limits();
  limits.max_document_bytes = 4;
  view.limits(limits);
  if(view.html("<p>too large</p>") ||
     view.last_error().find("exceeds byte limit") == std::string::npos)
  {
    return 1;
  }

  window.show();
  Fl::add_timeout(0.15, close_window, &window);
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
