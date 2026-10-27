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
#include "../FltkGedWindow.H"
#include "fltkcad/Fl_Cad_Terminal.H"
#include "bu/app.h"
#include "bu/log.h"
#include "bv/defines.h"
#include "dm/faceplate.h"
#include "ged.h"
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Menu_Bar.H>
#include <string>
#include <FL/Fl.H>
#include <FL/Fl_Preferences.H>
#include <cstdlib>
#include <iostream>

namespace {
void require(bool condition, const char* message)
{
  if (!condition) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

void menu_action(FltkGedWindow& window, const char* name)
{
  auto* menu = dynamic_cast<Fl_Menu_Bar*>(window.graphics_window().child(0));
  require(menu, "graphics menu missing");
  const auto* item = menu->find_item(name);
  require(item, "faceplate menu item missing");
  item->do_callback(menu);
}

Fl_Cad_View::InputEvent faceplate_event(Fl_Cad_View& view, double x, double y, int type)
{
  Fl_Cad_View::InputEvent event;
  event.type = type;
  event.button = FL_MIDDLE_MOUSE;
  event.x = static_cast<int>((x + 1) * view.w() * 0.5);
  event.y = static_cast<int>((1 - y) * view.h() * 0.5);
  return event;
}
}

int main(int, char** argv)
{
  bu_setprogname(argv[0]);
  if (!std::getenv("DISPLAY")) return 77;
  Fl_Preferences::file_access(Fl_Preferences::NONE);
  FltkGedWindow window(800, 300, Fl_Cad_View::Backend::OSMesa);
  window.show();
  Fl::check();
  require(!window.combined() && window.graphics_window().shown(),
          "default layout did not show separate command and graphics windows");
  require(window.view().window() == &window.graphics_window() &&
          window.terminal().window() == &window,
          "view and terminal were not assigned to their respective windows");
  require(window.view().w() == window.graphics_window().w(),
          "graphics view still reserves space for a sidebar");
  require(window.terminal().y() + window.terminal().h() == window.h(),
          "terminal still reserves a status bar for graphics information");
  auto& view = window.view();
  require(view.view()->gv_s->gv_view_params.draw && !view.faceplate_gui(),
          "view information or optional faceplate GUI defaults changed");
  require(view.view()->gv_s->gv_view_params.font_size == 0,
          "status font is not automatic by default");
  require(window.execute("view faceplate params font_size 12") == BRLCAD_OK &&
          view.view()->gv_s->gv_view_params.font_size == 12,
          "fixed status font command failed");
  require(window.execute("view faceplate params font_size invalid") == BRLCAD_ERROR &&
          view.view()->gv_s->gv_view_params.font_size == 12,
          "invalid status font changed the setting");
  window.terminal().resize(0, window.terminal().y(), 240, window.terminal().h());
  require(window.execute("view faceplate params font_size auto") == BRLCAD_OK &&
          view.view()->gv_s->gv_view_params.font_size == 0 &&
          ged_terminal_columns(window.gedp()) == static_cast<std::size_t>(window.terminal().columns()),
          "automatic status font or GED width synchronization failed");
  window.terminal().resize(0, window.terminal().y(), window.w(), window.terminal().h());
  const auto overlay = view.capture();
  require(overlay.has_value(), "could not capture view information overlay");
  menu_action(window, "&Modes/&Faceplate");
  const auto plain = view.capture();
  require(plain && plain->rgba != overlay->rgba,
          "view information was not drawn inside the graphics view");
  menu_action(window, "&Modes/Faceplate &GUI");
  const auto* gui = view.faceplate_gui();
  require(gui && gui->menu && gui->sliders && view.view()->gv_s->gv_view_params.draw,
          "faceplate GUI did not enable registered libdm controls");
  const auto controls = view.capture();
  require(controls && controls->rgba != overlay->rgba,
          "registered faceplate controls were not drawn");
  constexpr std::size_t zoom_in_row = 10;
  const double menu_y = gui->menu->top - gui->menu->row_height * (zoom_in_row + 0.5);
  const double scale = view.view()->gv_scale;
  require(view.process_input(faceplate_event(view, gui->menu->label_x, menu_y, FL_PUSH)) &&
          EQUAL(view.view()->gv_scale, scale), "menu probe changed the view");
  require(view.process_input(faceplate_event(view, gui->menu->label_x, menu_y, FL_RELEASE)) &&
          EQUAL(view.view()->gv_scale, scale * 0.5), "menu action did not zoom the view");
  point_t center_before, center_after;
  MAT_DELTAS_GET_NEG(center_before, view.view()->gv_center);
  const double slider_y = gui->sliders->top - gui->sliders->row_height * 0.5;
  view.process_input(faceplate_event(view, 0.5, slider_y, FL_PUSH));
  view.process_input(faceplate_event(view, 0.75, slider_y, FL_DRAG));
  view.process_input(faceplate_event(view, 0.75, slider_y, FL_RELEASE));
  MAT_DELTAS_GET_NEG(center_after, view.view()->gv_center);
  require(DIST_PNT_PNT_SQ(center_before, center_after) > SMALL_FASTF,
          "faceplate slider did not translate the view");
  menu_action(window, "&Modes/&Faceplate");
  require(!view.faceplate_gui() && !view.view()->gv_s->gv_view_params.draw,
          "hiding the faceplate left its controls active");
  menu_action(window, "&Modes/&Faceplate");
  require(view.faceplate_gui(), "showing the faceplate lost its GUI setting");
  menu_action(window, "&Modes/Faceplate &GUI");
  require(!view.faceplate_gui(), "faceplate GUI could not be disabled");
  window.terminal().input("pending");
  bu_log("terminal diagnostic hook test\n");
  Fl::wait(0.04);
  const bool logged = window.terminal().text().find("terminal diagnostic hook test") != std::string::npos;
  require(logged && window.terminal().input() == "pending",
          "library diagnostics did not reach the terminal or disturbed pending input");
  window.combined(true);
  Fl::check();
  require(!window.graphics_window().shown() && window.view().window() == &window,
          "combined mode did not reparent the graphics view");
  require(window.view().y() + window.view().h() == window.terminal().y(),
          "combined graphics and terminal panes do not meet");
  window.resize(window.x(), window.y(), 1000, 700);
  require(window.view().w() == window.w() && window.terminal().w() == window.w(),
          "combined panes did not track window resizing");
  window.combined(false);
  Fl::check();
  require(window.graphics_window().shown() &&
          window.view().window() == &window.graphics_window() &&
          window.terminal().h() == window.terminal().parent()->h(),
          "separating windows did not restore their panes");
  window.hide();
  require(!window.graphics_window().shown() && !window.shown(),
          "closing the command window left the graphics window visible");
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
