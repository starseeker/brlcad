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

#include "fltkcad/Fl_Ged_Accessibility.H"

#include <FL/Fl.H>
#include <FL/Fl_Browser_.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input_.H>
#include <FL/Fl_Valuator.H>

#include <string>

namespace {

bool has_text(const char* value)
{
  return value && *value;
}

void configure_group(Fl_Group& group)
{
  for (int index = 0; index < group.children(); ++index) {
    Fl_Widget* widget = group.child(index);
    if (!widget) {
      continue;
    }
    if (!has_text(widget->tooltip()) && has_text(widget->label())) {
      widget->copy_tooltip(widget->label());
    }
    if (auto* child = dynamic_cast<Fl_Group*>(widget)) {
      configure_group(*child);
    }
  }
}

void audit_group(const Fl_Group& group, fltkcad::AccessibilityAudit& audit)
{
  for (int index = 0; index < group.children(); ++index) {
    const Fl_Widget* widget = group.child(index);
    if (!widget) {
      continue;
    }
    ++audit.widgets_checked;
    const bool interactive = dynamic_cast<const Fl_Button*>(widget) ||
                             dynamic_cast<const Fl_Input_*>(widget) ||
                             dynamic_cast<const Fl_Browser_*>(widget) ||
                             dynamic_cast<const Fl_Choice*>(widget) ||
                             dynamic_cast<const Fl_Valuator*>(widget);
    if (interactive && !has_text(widget->label()) &&
        !has_text(widget->tooltip())) {
      audit.issues.emplace_back("Interactive widget has no label or tooltip");
    }
    if (widget->w() <= 0 || widget->h() <= 0) {
      audit.issues.emplace_back("Widget has no usable layout extent");
    }
    if (const auto* child = dynamic_cast<const Fl_Group*>(widget);
        child && !dynamic_cast<const Fl_Browser_*>(widget)) {
      audit_group(*child, audit);
    }
  }
}

} // namespace

void fltkcad::configure_accessibility(Fl_Group& root)
{
  Fl::option(Fl::OPTION_VISIBLE_FOCUS, true);
  configure_group(root);
}

fltkcad::AccessibilityAudit fltkcad::audit_accessibility(const Fl_Group& root)
{
  AccessibilityAudit audit;
  audit_group(root, audit);
  return audit;
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
