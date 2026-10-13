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

#include "fltkcad/CommandHistory.H"

#include <cstdlib>
#include <iostream>
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

} // namespace

int
main(int, char** argv)
{
  bu_setprogname(argv[0]);
  fltkcad::CommandHistory history(3);
  history.add("one");
  history.add("one");
  history.add("two");
  history.add("three");
  history.add("four");

  require(history.size() == 3, "history limit or duplicate suppression failed");
  require(history.at(0) == "two", "old history was not evicted");
  require(history.previous("draft") == "four", "history previous failed");
  require(history.previous("ignored") == "three", "history traversal failed");
  require(history.next() == "four", "history next failed");
  require(history.next() == "draft", "draft restoration failed");
  require(!history.next(), "history advanced beyond draft");

  require(history.consolidate(0, 2), "valid consolidation failed");
  require(history.size() == 2, "consolidation size is wrong");
  require(history.at(0) == "two three", "consolidation contents are wrong");
  require(!history.consolidate(1, 1), "empty consolidation was accepted");
  require(history.consolidate(0, 2, "set value {first\nsecond}") &&
          history.at(0) == "set value {first\nsecond}",
          "consolidation changed multiline Tcl history");

  history.limit(1);
  require(history.limit() == 1 && history.size() == 1,
          "reducing the history limit did not evict old commands");

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
