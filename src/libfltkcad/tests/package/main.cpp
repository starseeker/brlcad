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
#include "dm/faceplate.h"

#include <fltkcad/CommandHistory.H>
#include <fltkcad/Fl_Cad_Quad_View.H>
#include <fltkcad/Fl_Cad_Terminal.H>
#include <fltkcad/Fl_Cad_View.H>
#include <fltkcad/Fl_Ged_Accessibility.H>
#include <fltkcad/Fl_Ged_Attributes.H>
#include <fltkcad/Fl_Ged_Combination_Editor.H>
#include <fltkcad/Fl_Ged_Create_Primitive.H>
#include <fltkcad/Fl_Ged_Help.H>
#include <fltkcad/Fl_Ged_Hierarchy.H>
#include <fltkcad/Fl_Ged_Interaction.H>
#include <fltkcad/Fl_Ged_Preferences.H>
#include <fltkcad/Fl_Ged_Primitive_Editor.H>
#include <fltkcad/Fl_Ged_Render_Tools.H>
#include <fltkcad/Fl_Ged_Search.H>
#include <fltkcad/Fl_Ged_Sketch_Editor.H>
#include <fltkcad/Fl_Ged_Specialist_Tools.H>
#include <fltkcad/Fl_Ged_Tcl.H>
#include <fltkcad/Fl_Ged_Transform_Editor.H>
#include <fltkcad/Fl_Ged_View_Controls.H>
#include <fltkcad/Fl_Litehtml_View.H>

#include <type_traits>

static_assert(std::is_base_of_v<Fl_Group, Fl_Litehtml_View>);

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  struct dm_faceplate_controls menu;
  dm_faceplate_controls_init(&menu, DM_FACEPLATE_MENU);
  fltkcad::CommandHistory history(1);
  history.add("draw all");
  return history.size() == 1 && history.at(0) == "draw all" ? 0 : 1;
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
