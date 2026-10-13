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

#ifndef FLTKCAD_FL_GED_ACCESSIBILITY_H
#define FLTKCAD_FL_GED_ACCESSIBILITY_H

#include "fltkcad/defines.h"

#include <cstddef>
#include <string>
#include <vector>

class Fl_Group;

namespace fltkcad {

struct FLTKCAD_EXPORT AccessibilityAudit
{
  std::size_t widgets_checked = 0;
  std::vector<std::string> issues;
};

FLTKCAD_EXPORT void configure_accessibility(Fl_Group& root);
FLTKCAD_EXPORT AccessibilityAudit audit_accessibility(const Fl_Group& root);

} // namespace fltkcad

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
