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

#include "Uri.h"

#include <filesystem>
#include <string>

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  using fltkcad::detail::file_uri;
  using fltkcad::detail::local_path;
  using fltkcad::detail::resolve_uri;

  if(resolve_uri("style.css", "memory://manual/index.html") !=
       "memory://manual/style.css" ||
     resolve_uri("#target", "memory://manual/index.html?mode=1") !=
       "memory://manual/index.html#target" ||
     resolve_uri("?mode=2", "memory://manual/index.html?mode=1") !=
       "memory://manual/index.html?mode=2" ||
     resolve_uri("//cdn.example/image.png",
                 "https://example/manual/index.html") !=
       "https://cdn.example/image.png" ||
     local_path("https://example/manual/index.html"))
  {
    return 1;
  }

  const auto path = std::filesystem::current_path() / "manual" / "index.html";
  const auto round_trip = local_path(file_uri(path));
  if(!round_trip || round_trip->lexically_normal() != path.lexically_normal())
  {
    return 1;
  }
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
