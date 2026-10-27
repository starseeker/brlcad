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
#include "fltkcad/Fl_Ged_Combination_Editor.H"
#include "bu/app.h"
#include "bu/file.h"
#include "bv/util.h"
#include "ged.h"
#include "ged/edit.h"
#include "dm/faceplate.h"
#include <FL/Fl.H>
#include <FL/Fl_Preferences.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Menu_Bar.H>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace {
constexpr double kViewTolerance = 1.0e-8;
constexpr double kAngleTolerance = 1.0e-6;
constexpr point_t kGeometryCenter = {10000, 20000, -30000};
constexpr const char* kGeometry = "fltkged_application_primary.s";
constexpr const char* kOtherGeometry = "fltkged_application_secondary.s";

void require(bool condition, const std::string& message)
{
  if (!condition) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

class DatabaseFixture
{
public:
  explicit DatabaseFixture(const char* source)
  {
    char filename[MAXPATHLEN] = {};
    FILE* file = bu_temp_file(filename, sizeof(filename));
    require(file, "could not create temporary database fixture");
    require(std::fclose(file) == 0, "could not close temporary database fixture");
    path_ = filename;
    std::filesystem::copy_file(source, path_, std::filesystem::copy_options::overwrite_existing);
  }
  ~DatabaseFixture() { bu_file_delete(path_.c_str()); }
  const std::string& path() const { return path_; }
private:
  std::string path_;
};

void execute(FltkGedWindow& window, const std::string& command)
{
  require(window.execute(command) == BRLCAD_OK, "command failed: " + command);
}

void check_clear(FltkGedWindow& window, const std::string& command)
{
  auto& terminal = window.terminal();
  terminal.append("output to clear\nmore output\n");
  terminal.input(command);
  terminal.submit();
  const std::string screen = terminal.text();
  require(screen == terminal.prompt() && terminal.input().empty() &&
          terminal.output().insert_position() == static_cast<int>(terminal.prompt().size()),
          "clear did not leave a clean screen and usable prompt for '" + command +
          "': " + screen);
  require(terminal.history_at(terminal.history_size() - 1) == command,
          "clear erased command history");
}

template <typename Widget>
Widget* find_widget(Fl_Group& group, const char* label)
{
  for (int index = 0; index < group.children(); ++index) {
    Fl_Widget* child = group.child(index);
    if (auto* match = dynamic_cast<Widget*>(child)) {
      if (match->label() && std::string(match->label()) == label) return match;
    }
    if (auto* nested = dynamic_cast<Fl_Group*>(child)) {
      if (auto* match = find_widget<Widget>(*nested, label)) return match;
    }
  }
  return nullptr;
}

void cancel_quit(void*)
{
  require(Fl::modal(), "quit did not prompt for unsaved changes");
  Fl::modal()->do_callback();
}

void check_editor_guards(FltkGedWindow& window, const std::string& path)
{
  execute(window, std::string("comb fltkged_application.c u ") + kGeometry);
  execute(window, "tcl fltkged selection fltkged_application.c");
  auto* menu = dynamic_cast<Fl_Menu_Bar*>(window.child(0));
  require(menu, "application menu missing");
  const auto* item = menu->find_item("&Edit/Edit Selected &Combination...");
  require(item, "combination editor menu missing");
  item->do_callback(menu);
  Fl_Ged_Combination_Editor* editor = nullptr;
  for (Fl_Window* current = Fl::first_window(); current; current = Fl::next_window(current)) {
    if ((editor = dynamic_cast<Fl_Ged_Combination_Editor*>(current))) break;
  }
  require(editor, "combination editor did not open");
  auto* shader = find_widget<Fl_Input>(*editor, "Shader");
  require(shader, "combination shader input missing");
  shader->value("plastic");
  shader->do_callback();
  require(editor->dirty(), "combination edit was not marked unsaved");
  auto* original_database = window.gedp()->dbip;
  execute(window, "opendb");
  require(editor->dirty(), "querying the database discarded unsaved edits");
  require(window.execute("closedb") == BRLCAD_ERROR &&
          window.gedp()->dbip == original_database && editor->dirty(),
          "closedb bypassed the unsaved editor guard");
  require(window.execute("opendb \"" + path + "\"") == BRLCAD_ERROR &&
          window.gedp()->dbip == original_database && editor->dirty(),
          "opendb bypassed the unsaved editor guard");
  Fl::add_timeout(0.01, cancel_quit);
  execute(window, "quit");
  require(window.shown() && window.graphics_window().shown() && editor->dirty(),
          "cancelling quit closed the application or discarded edits");
  require(editor->save(), "could not save the edited combination");
}

void check_edit_session(FltkGedWindow& window)
{
  execute(window, std::string("sed ") + kGeometry);
  execute(window, "tra -m 1 0 0");
  auto* original_database = window.gedp()->dbip;
  require(window.execute("closedb") == BRLCAD_ERROR &&
          window.gedp()->dbip == original_database,
          "closedb discarded a command-driven edit");
  Fl::add_timeout(0.01, cancel_quit);
  execute(window, "quit");
  require(window.shown() && ged_edit_session_mode(window.gedp()) == GED_EDIT_PRIMITIVE,
          "cancelling quit discarded a command-driven edit");

  auto* menu = dynamic_cast<Fl_Menu_Bar*>(window.child(0));
  require(menu, "application menu missing");
  const auto* toggle = menu->find_item("&Modes/Faceplate &GUI");
  require(toggle, "faceplate GUI menu missing");
  toggle->do_callback(menu);
  const auto* gui = window.view().faceplate_gui();
  require(gui && gui->menu, "faceplate GUI did not enable");
  struct dm_faceplate_item action{};
  require(gui->menu->item(gui->menu->data, 1, &action) && action.label &&
          std::string(action.label) == "Accept", "faceplate did not expose the active edit actions");
  gui->menu->select(gui->menu->data, 4, 0, 1);
  require(ged_edit_session_mode(window.gedp()) == GED_EDIT_VIEW,
          "faceplate Reject did not end the shared edit");
  toggle->do_callback(menu);
}
}

int main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  if (!std::getenv("DISPLAY")) return 77;
  require(argc == 2, "a database fixture is required");
  Fl_Preferences::file_access(Fl_Preferences::NONE);
  DatabaseFixture database(argv[1]);
  FltkGedWindow window(800, 300, Fl_Cad_View::Backend::OSMesa);
  window.show();
  Fl::check();
  auto* view = window.view().view();
  vect_t aet;
  bv_view_get_aet(view, aet);
  require(NEAR_EQUAL(aet[0], 35, kAngleTolerance) &&
          NEAR_EQUAL(aet[1], 25, kAngleTolerance) && NEAR_ZERO(aet[2], kAngleTolerance),
          "startup view does not match MGED's 35/25 orientation: " +
          std::to_string(aet[0]) + ", " + std::to_string(aet[1]) + ", " +
          std::to_string(aet[2]));
  check_clear(window, "clear");
  require(window.open_database(database.path()), "database fixture did not open");
  check_clear(window, "  clear  ");
  if (window.command_language() == FltkGedWindow::CommandLanguage::Tcl) {
    const std::size_t history_begin = window.terminal().history_size();
    window.terminal().input("set history_value {first");
    window.terminal().submit();
    window.terminal().input("second}");
    window.terminal().submit();
    require(window.terminal().history_size() == history_begin + 1 &&
            window.terminal().history_at(history_begin) == "set history_value {first\nsecond}",
            "multiline Tcl history did not preserve its script");
    require(window.execute_prompt("set prompt_value {first") == GED_MORE,
            "Tcl prompt did not accept an incomplete command");
    require(window.execute_prompt("second}; puts $prompt_value") == BRLCAD_OK,
            "Tcl prompt did not complete a multiline command");
    require(window.execute_prompt("in prompt.s sph") == GED_MORE,
            "Tcl prompt did not preserve GED_MORE");
    require(window.execute_prompt("0 0 0 1") == BRLCAD_OK,
            "Tcl prompt did not finish GED_MORE input");
    require(window.execute_prompt("kill prompt.s; expr {[llength [list a b]] == 2}") == BRLCAD_OK,
            "Tcl substitutions did not evaluate in the command prompt");
    require(window.execute_prompt("error {expected prompt error}") == BRLCAD_ERROR,
            "Tcl errors were not reported");
    window.command_language(FltkGedWindow::CommandLanguage::Ged);
    require(window.execute_prompt("size 1000") == BRLCAD_OK,
            "raw GED prompt did not execute");
    window.command_language(FltkGedWindow::CommandLanguage::Tcl);
    require(window.execute_prompt("expr {$prompt_value eq {first\nsecond}}") == BRLCAD_OK,
            "switching shells lost Tcl state");
  }
  execute(window, std::string("in ") + kGeometry + " sph 10000 20000 -30000 125");
  execute(window, std::string("in ") + kOtherGeometry + " sph 50000 20000 -30000 250");
  const double initial_scale = bv_view_get_scale(view);
  execute(window, std::string("e ") + kGeometry);
  point_t center;
  bv_view_get_center_vec(view, center);
  const double fitted_scale = bv_view_get_scale(view);
  require(DIST_PNT_PNT_SQ(center, kGeometryCenter) < kViewTolerance &&
          !NEAR_EQUAL(fitted_scale, initial_scale, kViewTolerance),
          "first draw did not automatically frame the offset geometry: " +
          std::to_string(center[0]) + ", " + std::to_string(center[1]) + ", " +
          std::to_string(center[2]) + "; scale " + std::to_string(fitted_scale));
  execute(window, "autoview");
  require(NEAR_EQUAL(bv_view_get_scale(view), fitted_scale, kViewTolerance),
          "automatic framing differs from the autoview command");
  execute(window, "size 2000");
  execute(window, std::string("draw ") + kOtherGeometry);
  require(NEAR_EQUAL(bv_view_get_size(view), 2000, kViewTolerance),
          "additional draw reset the user's view");
  check_clear(window, "clear");
  require(ged_who_argc(window.gedp()) == 2, "clear erased displayed geometry");
  execute(window, "Z");
  require(BU_PTBL_LEN(bv_view_objs(view, BV_DB_OBJS)) == 0,
          "Z left recycled geometry in the view's object table");
  execute(window, std::string("draw -R ") + kGeometry);
  require(NEAR_EQUAL(bv_view_get_size(view), 2000, kViewTolerance),
          "draw -R did not preserve the view on an empty display");
  execute(window, std::string("B ") + kOtherGeometry);
  bv_view_get_center_vec(view, center);
  require(NEAR_EQUAL(center[0], 50000, kViewTolerance),
          "B did not frame the replacement geometry");
  execute(window, "Z");
  execute(window, "size 2000");
  execute(window, std::string("draw -m0R ") + kGeometry);
  require(NEAR_EQUAL(bv_view_get_scale(view), fitted_scale, kViewTolerance),
          "an option argument containing R was treated as a resize flag");
  const bool tcl_available = window.execute("tcl expr 1") == BRLCAD_OK;
  if (tcl_available) {
    check_clear(window, "tcl ged clear");
    execute(window, "Z");
    execute(window, std::string("tcl ged draw ") + kGeometry);
    bv_view_get_center_vec(view, center);
    require(DIST_PNT_PNT_SQ(center, kGeometryCenter) < kViewTolerance &&
            NEAR_EQUAL(bv_view_get_scale(view), fitted_scale, kViewTolerance),
            "Tcl drawing did not apply automatic framing");
  }
  const double scale_before_error = bv_view_get_scale(view);
  require(window.execute(std::string("draw -L -1 ") + kGeometry) != BRLCAD_OK &&
          NEAR_EQUAL(bv_view_get_scale(view), scale_before_error, kViewTolerance),
          "failed draw changed the camera");
  execute(window, std::string("search . -name ") + kGeometry +
                  " -exec echo frontend_search_exec {} ;");
  const bool executed = window.terminal().text().find("frontend_search_exec") != std::string::npos;
  require(executed, "search -exec did not execute its command");
  if (tcl_available) check_editor_guards(window, database.path());
  check_edit_session(window);
  execute(window, "close");
  require(!window.gedp()->dbip && std::string(window.label()).find(database.path()) == std::string::npos,
          "close did not update the application database state");
  execute(window, "open \"" + database.path() + "\"");
  require(window.gedp()->dbip && std::string(window.label()).find(database.path()) != std::string::npos,
          "open did not update the application database state");
  for (const char* command : {"quit", "q", "exit"}) {
    window.show();
    require(window.execute(std::string(command) + " unexpected") == GED_HELP && window.shown(),
            "invalid quit arguments closed the application");
    execute(window, command);
    require(!window.shown() && !window.graphics_window().shown(),
            std::string(command) + " left application windows open");
  }
  if (tcl_available) {
    window.show();
    execute(window, "tcl ged quit");
    require(!window.shown(), "Tcl quit did not close the application");
  }
  window.show();
  require(window.open_database(database.path()), "database replacement failed");
  check_clear(window, "clear");
  execute(window, "quit");
  require(!window.shown(), "database replacement lost the quit callback");
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
