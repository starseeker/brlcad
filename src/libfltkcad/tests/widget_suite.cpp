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
#include "vmath.h"

#include "fltkcad/Fl_Ged_Combination_Editor.H"
#include "fltkcad/Fl_Ged_Accessibility.H"
#include "fltkcad/Fl_Ged_Create_Primitive.H"
#include "fltkcad/Fl_Ged_Primitive_Editor.H"
#include "fltkcad/Fl_Ged_Preferences.H"
#include "fltkcad/Fl_Ged_Interaction.H"
#include "fltkcad/Fl_Ged_Render_Tools.H"
#include "fltkcad/Fl_Ged_Search.H"
#include "fltkcad/Fl_Ged_Sketch_Editor.H"
#include "fltkcad/Fl_Ged_Specialist_Tools.H"
#include "fltkcad/Fl_Ged_Tcl.H"
#include "fltkcad/Fl_Ged_Transform_Editor.H"
#include "fltkcad/Fl_Ged_View_Controls.H"
#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "dm.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/view.h"
#include "bv/util.h"
#include "bu/str.h"
#include "raytrace.h"
#include "rt/search.h"
}

namespace {

constexpr double kViewTolerance = 1.0e-8;

bool delete_fixture_file(const std::filesystem::path& path)
{
  const auto filename = path.string();
  if (!bu_file_exists(filename.c_str(), nullptr) ||
      bu_file_delete(filename.c_str())) {
    return true;
  }
  std::cerr << "Unable to delete test fixture: " << filename << '\n';
  return false;
}

class Fixture
{
public:
  explicit Fixture(const std::filesystem::path& source)
  {
    const auto suffix = std::chrono::high_resolution_clock::now()
                          .time_since_epoch()
                          .count();
    path_ = std::filesystem::temp_directory_path() /
            ("fltkged-widgets-" + std::to_string(suffix) + ".g");
    std::filesystem::copy_file(source, path_);
    ged_ = ged_open("db", path_.string().c_str(), 1);
  }

  ~Fixture()
  {
    if (ged_) {
      ged_close(ged_);
    }
    delete_fixture_file(path_);
  }

  Fixture(const Fixture&) = delete;
  Fixture& operator=(const Fixture&) = delete;

  struct ged* get() const noexcept { return ged_; }
  const std::filesystem::path& path() const noexcept { return path_; }

private:
  std::filesystem::path path_;
  struct ged* ged_ = nullptr;
};

std::string first_object(struct ged* gedp, int flags, int excluded_type = -1)
{
  struct directory** entries = nullptr;
  const std::size_t count = db_ls(gedp->dbip, flags, nullptr, &entries);
  std::string result;
  for (std::size_t index = 0; index < count; ++index) {
    if (entries[index] && entries[index]->d_namep &&
        entries[index]->d_minor_type != excluded_type) {
      result = entries[index]->d_namep;
      break;
    }
  }
  if (entries) {
    bu_free(entries, "widget test entries");
  }
  return result;
}

int fail(const std::string& message)
{
  std::cerr << message << '\n';
  return 1;
}

template <typename Widget>
Widget* find_widget(Fl_Group& parent, const std::string& label)
{
  for (int index = 0; index < parent.children(); ++index) {
    Fl_Widget* child = parent.child(index);
    if (auto* match = dynamic_cast<Widget*>(child);
        match && match->label() && label == match->label()) {
      return match;
    }
    if (auto* group = dynamic_cast<Fl_Group*>(child)) {
      if (auto* match = find_widget<Widget>(*group, label)) {
        return match;
      }
    }
  }
  return nullptr;
}

bool bump_numeric_input(Fl_Group& parent)
{
  for (int index = 0; index < parent.children(); ++index) {
    Fl_Widget* child = parent.child(index);
    if (auto* input = dynamic_cast<Fl_Input*>(child)) {
      std::istringstream text(input->value() ? input->value() : "");
      std::vector<double> values;
      double value = 0.0;
      while (text >> value) {
        values.push_back(value);
      }
      if (!values.empty() && text.eof()) {
        if (EQUAL(std::trunc(values.front()), values.front())) {
          values.front() += 1.0;
        } else {
          values.front() +=
            std::max(0.001, std::abs(values.front()) * 0.01);
        }
        std::ostringstream replacement;
        replacement << std::setprecision(12);
        for (std::size_t value_index = 0; value_index < values.size();
             ++value_index) {
          if (value_index) {
            replacement << ' ';
          }
          replacement << values[value_index];
        }
        input->value(replacement.str().c_str());
        return true;
      }
    }
    if (auto* group = dynamic_cast<Fl_Group*>(child);
        group && bump_numeric_input(*group)) {
      return true;
    }
  }
  return false;
}

} // namespace

int main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  if (argc != 3) {
    return fail("usage: widget_suite BACKEND DATABASE");
  }
  const Fl_Cad_View::Backend backend =
    std::string(argv[1]) == "opengl" ? Fl_Cad_View::Backend::SystemOpenGL
                                     : Fl_Cad_View::Backend::OSMesa;

  Fixture fixture(argv[2]);
  struct ged* gedp = fixture.get();
  if (!gedp || !gedp->dbip) {
    return fail("unable to open test database copy");
  }

  const char* make_sketch[] = {"make", "fltkged_test_sketch", "sketch"};
  bu_vls_trunc(gedp->ged_result_str, 0);

  const char* make_bot[] = {"make", "fltkged_test_bot", "bot"};
  if (ged_exec(gedp, 3, make_bot) != BRLCAD_OK) {
    return fail(std::string("unable to create test BOT: ") +
                bu_vls_cstr(gedp->ged_result_str));
  }
  bu_vls_trunc(gedp->ged_result_str, 0);

  const char* make_nmg[] = {"make", "fltkged_test_nmg", "nmg"};
  if (ged_exec(gedp, 3, make_nmg) != BRLCAD_OK) {
    return fail(std::string("unable to create test NMG: ") +
                bu_vls_cstr(gedp->ged_result_str));
  }
  bu_vls_trunc(gedp->ged_result_str, 0);
  if (ged_exec(gedp, 3, make_sketch) != BRLCAD_OK) {
    return fail(std::string("unable to create test sketch: ") +
                bu_vls_cstr(gedp->ged_result_str));
  }
  bu_vls_trunc(gedp->ged_result_str, 0);

  {
    Fl_Ged_Create_Primitive create;
    create.ged(gedp);
    auto* name = find_widget<Fl_Input>(create, "Name");
    auto* type = find_widget<Fl_Choice>(create, "Type");
    if (!name || !type) {
      return fail("create primitive controls were not found");
    }
    name->value("fltkged_test_sphere");
    int sphere_index = -1;
    for (int index = 0; index < type->size(); ++index) {
      if (type->text(index) && std::string(type->text(index)) == "sph") {
        sphere_index = index;
        break;
      }
    }
    type->value(sphere_index);
    if (sphere_index < 0 || !create.create() ||
        !db_lookup(gedp->dbip, "fltkged_test_sphere", LOOKUP_QUIET)) {
      return fail("create primitive dialog failed: " + create.last_error());
    }
  }

  const std::string primitive = "fltkged_test_sphere";
  const std::string combination = first_object(gedp, DB_LS_COMB);
  if (primitive.empty() || combination.empty()) {
    return fail("test database does not contain required object types");
  }

  {
    Fl_Ged_Primitive_Editor editor(backend);
    if (editor.open(gedp, "missing-object") || editor.last_error().empty()) {
      return fail("primitive editor did not reject a missing object");
    }
    if (!editor.open(gedp, primitive)) {
      return fail("primitive editor failed to open " + primitive + ": " +
                  editor.last_error());
    }
    editor.show();
    Fl::check();
    auto* operations = find_widget<Fl_Choice>(editor, "Operation");
    bool applied = false;
    for (int index = 0; operations && index < operations->size(); ++index) {
      operations->value(index);
      operations->do_callback();
      bump_numeric_input(editor);
      if (editor.apply()) {
        applied = true;
        break;
      }
    }
    if (!applied || !editor.dirty() || !editor.save() || editor.dirty()) {
      return fail("primitive editor save/reload cycle failed: " +
                  editor.last_error());
    }
    editor.hide();
  }

  {
    Fl_Ged_Combination_Editor editor;
    if (!editor.open(gedp, combination)) {
      return fail("combination editor failed to open " + combination + ": " +
                  editor.last_error());
    }
    if (editor.dirty()) {
      return fail("combination editor is dirty immediately after loading");
    }
    auto* shader = find_widget<Fl_Input>(editor, "Shader");
    auto* add_union = find_widget<Fl_Button>(editor, "Add Union");
    auto* add_subtract = find_widget<Fl_Button>(editor, "Add Subtract");
    auto* add_intersect = find_widget<Fl_Button>(editor, "Add Intersect");
    if (!shader || !add_union || !add_subtract || !add_intersect) {
      return fail("combination editor member controls were not found");
    }
    shader->do_callback();
    if (!editor.dirty() || !editor.save() || editor.dirty()) {
      return fail("combination editor save cycle failed: " +
                  editor.last_error());
    }
  }

  {
    Fl_Ged_Sketch_Editor editor(backend);
    if (!editor.open(gedp, "fltkged_test_sketch")) {
      return fail("sketch editor failed to open its fixture: " +
                  editor.last_error());
    }
    editor.show();
    Fl::check();
    auto* u = find_widget<Fl_Input>(editor, "U");
    auto* v = find_widget<Fl_Input>(editor, "V");
    auto* add = find_widget<Fl_Button>(editor, "Add Vertex");
    auto* add_arc = find_widget<Fl_Button>(editor, "Add Arc...");
    auto* add_bezier = find_widget<Fl_Button>(editor, "Add Bezier...");
    auto* add_nurb = find_widget<Fl_Button>(editor, "Add NURB...");
    auto* start = find_widget<Fl_Input>(editor, "Start");
    auto* end = find_widget<Fl_Input>(editor, "End");
    auto* add_line = find_widget<Fl_Button>(editor, "Add Line");
    if (!u || !v || !add || !add_arc || !add_bezier || !add_nurb || !start ||
        !end || !add_line) {
      return fail("sketch editor vertex controls were not found");
    }
    u->value("0");
    v->value("0");
    add->do_callback();
    u->value("10");
    v->value("3");
    add->do_callback();
    u->value("20");
    v->value("5");
    add->do_callback();
    if (!editor.constrain_vertices(
          Fl_Ged_Sketch_Editor::VertexConstraint::Horizontal, 0, 1) ||
        !editor.move_vertices({0, 1}, 1.0, 1.0)) {
      return fail("sketch constraint or vertex-list move failed: " +
                  editor.last_error());
    }
    start->value("0");
    end->value("1");
    add_line->do_callback();
    if (!editor.split_segment(0, 0.5) || !editor.reverse_segment(0)) {
      return fail("advanced sketch segment edit failed: " +
                  editor.last_error());
    }
    if (editor.constrain_vertices(
          Fl_Ged_Sketch_Editor::VertexConstraint::Distance, 0, 99, 1.0) ||
        editor.last_error().empty()) {
      return fail("sketch editor accepted an invalid constraint");
    }
    if (!editor.dirty() || !editor.save() || editor.dirty()) {
      return fail("sketch editor save/reload cycle failed: " +
                  editor.last_error());
    }
    editor.hide();
  }

  {
    Fl_Ged_Search search;
    search.ged(gedp);
    if (!search.search()) {
      return fail("database search failed: " + search.last_error());
    }
  }

  {
    Fl_Ged_Tcl interpreter;
    interpreter.ged(gedp);
    if (interpreter.available()) {
      if (interpreter.evaluate("set value 7; expr {$value + 2}") != 0 ||
          interpreter.result() != "9") {
        return fail("embedded Tcl evaluation failed: " +
                    interpreter.last_error());
      }
      if (interpreter.evaluate("ged ls") != 0 ||
          interpreter.result().find(primitive) == std::string::npos) {
        return fail("embedded Tcl GED dispatcher failed: " +
                    interpreter.last_error());
      }
      if (interpreter.evaluate("expr {[string first {Usage: status} [help status]] >= 0 && "
                               "[lsearch -exact [?] sed] >= 0 && "
                               "[lsearch -exact [::ged::commands] sed] >= 0}") != 0 ||
          interpreter.result() != "1") {
        return fail("classic help and installed command inventory disagree: " + interpreter.last_error());
      }
      struct ged* unit_context = ged_open("inmem", "0", 0);
      if (!unit_context) return fail("could not create unit trace fixture");
      const char* inches[] = {"units", "in"};
      const bool units_set = ged_exec(unit_context, 2, inches) == BRLCAD_OK;
      interpreter.ged(unit_context);
      const bool units_match = units_set && interpreter.evaluate(
        "expr {abs($local2base - 25.4) < 1e-12 && abs($base2local * $local2base - 1) < 1e-12}") == 0 &&
        interpreter.result() == "1";
      interpreter.ged(nullptr);
      const bool closed_units_fail = interpreter.evaluate("set local2base") != 0;
      interpreter.ged(gedp);
      ged_close(unit_context);
      if (!units_match || !closed_units_fail) return fail("Tcl unit variables did not follow database binding");
      interpreter.ui_handler(
        [](const std::vector<std::string>& arguments, std::string& result) {
          result = arguments.empty() ? "none" : arguments.front();
          return 0;
        });
      if (interpreter.evaluate("fltkged accessibility") != 0 ||
          interpreter.result() != "accessibility") {
        return fail("embedded Tcl FLTK bridge failed: " +
                    interpreter.last_error());
      }
    }
    if (interpreter.source("/a/path/that/does/not/exist.tcl") == 0 ||
        interpreter.last_error().empty()) {
      return fail("Tcl source did not reject a missing script");
    }
  }

  {
    Fl_Ged_Specialist_Tools tools;
    tools.ged(gedp);
    std::string command;
    tools.command_handler([&command](const std::string& value) {
      command = value;
      return BRLCAD_OK;
    });
    tools.selection("fltkged_test_bot");
    if (!tools.inspect_bot() ||
        command != "bot check solid \"fltkged_test_bot\"") {
      return fail("BOT inspection panel did not issue its validated command: " +
                  tools.last_error());
    }
    if (!tools.repair_bot() ||
        command != "bot repair -o \"fltkged_test_bot.repaired\" \"fltkged_test_bot\"") {
      return fail("BOT repair panel did not issue its validated command: " +
                  tools.last_error());
    }
    tools.selection("fltkged_test_nmg");
    if (!tools.fix_nmg_normals() ||
        command != "nmg_fix_normals \"fltkged_test_nmg\"") {
      return fail("NMG repair panel did not issue its validated command: " +
                  tools.last_error());
    }
    tools.selection(combination);
    auto* shader = find_widget<Fl_Input>(tools, "Shader");
    auto* color = find_widget<Fl_Input>(tools, "RGB");
    auto* inherit = find_widget<Fl_Input>(tools, "Inherit");
    if (!shader || !color || !inherit) {
      return fail("material panel controls were not found");
    }
    shader->value("plastic");
    color->value("10 20 30");
    inherit->value("n");
    if (!tools.apply_material() || command.find("mater -s") != 0) {
      return fail("material panel did not issue its validated command: " +
                  tools.last_error());
    }
    std::filesystem::path density_path = fixture.path();
    density_path += ".density";
    std::filesystem::path export_path = fixture.path();
    export_path += ".export.density";
    if (!delete_fixture_file(density_path) || !delete_fixture_file(export_path)) {
      return fail("unable to clear density test fixtures");
    }
    {
      std::ofstream density_file(density_path);
      density_file << "1 1.0 test material\n";
      if (!density_file) {
        return fail("unable to create the density test fixture");
      }
    }
    if (!tools.validate_densities(density_path.string()) ||
        command.find("mater -d validate") != 0 ||
        !tools.import_densities(density_path.string()) ||
        command.find("mater -d import -v") != 0 ||
        !tools.export_densities(export_path.string()) ||
        command.find("mater -d export") != 0 || !tools.clear_densities() ||
        command != "mater -d clear") {
      delete_fixture_file(density_path);
      return fail("density tools did not issue validated commands: " +
                  tools.last_error());
    }
    if (!delete_fixture_file(density_path)) {
      return fail("unable to delete the density test fixture");
    }
    if (tools.validate_densities(density_path.string()) ||
        tools.last_error().empty()) {
      return fail("density tools accepted a missing input file");
    }
    tools.selection(primitive);
    if (tools.apply_material() || tools.last_error().empty()) {
      return fail("material panel accepted a non-combination object");
    }
    fltkcad::configure_accessibility(tools);
    const auto audit = fltkcad::audit_accessibility(tools);
    if (audit.widgets_checked == 0 || !audit.issues.empty()) {
      return fail("specialist panel accessibility audit failed");
    }
  }

  {
    Fl_Ged_Transform_Editor editor;
    editor.ged(gedp);
    editor.command_handler([gedp](const std::string& text) {
      std::vector<char> input(text.begin(), text.end());
      input.push_back('\0');
      const std::size_t max_args = 32;
      char* command_args[max_args];
      int command_count = bu_argv_from_string(command_args, max_args, input.data());
      return ged_edit_session_exec(gedp, command_count, const_cast<const char**>(command_args));
    });
    editor.selection("/" + primitive);
    auto* scale = find_widget<Fl_Input>(editor, "Uniform scale");
    auto* apply_scale = find_widget<Fl_Button>(editor, "Scale");
    auto* accept = find_widget<Fl_Button>(editor, "Accept");
    auto* reject = find_widget<Fl_Button>(editor, "Reject");
    if (!scale || !apply_scale || !accept || !reject)
      return fail("transform editor controls were not found");
    scale->value("1.25");
    apply_scale->do_callback();
    if (ged_edit_session_mode(gedp) != GED_EDIT_PRIMITIVE)
      return fail("transform editor did not enter primitive editing");
    reject->do_callback();
    if (ged_edit_session_mode(gedp) != GED_EDIT_VIEW)
      return fail("transform editor did not reject its working edit");

    const char* comb[] = {"comb", "fltkged_test_transform.c", "u", primitive.c_str()};
    if (ged_exec(gedp, 4, comb) != BRLCAD_OK)
      return fail("unable to create transform fixture");
    editor.selection("/fltkged_test_transform.c/" + primitive);
    auto* translation = find_widget<Fl_Input>(editor, "Translation");
    auto* translate = find_widget<Fl_Button>(editor, "Translate");
    translation->value("1 2 3");
    translate->do_callback();
    struct rt_edit* working = ged_edit_session_edit(gedp);
    if (ged_edit_session_mode(gedp) != GED_EDIT_INSTANCE || !working ||
        !NEAR_EQUAL(working->model_changes[3], gedp->dbip->dbi_local2base, kViewTolerance))
      return fail("transform editor did not preview the selected member translation");
    accept->do_callback();
    if (ged_edit_session_mode(gedp) != GED_EDIT_VIEW)
      return fail("transform editor did not accept its working edit");
  }

  {
    Fl_Ged_Render_Tools tools;
    std::string command;
    tools.command_handler([&command](const std::string& value) {
      command = value;
      return BRLCAD_OK;
    });
    auto* render = find_widget<Fl_Button>(tools, "Render in View");
    if (!render) {
      return fail("render tool controls were not found");
    }
    render->do_callback();
    if (command.find("ert") != 0) {
      return fail("render tools did not issue an embedded raytrace command");
    }
  }

  {
    Fl_Ged_Preferences preferences;
    Fl_Ged_Preferences::Settings changed;
    bool applied = false;
    preferences.changed_handler(
      [&changed, &applied](const Fl_Ged_Preferences::Settings& value) {
        changed = value;
        applied = true;
      });
    auto* first = find_widget<Fl_Input>(preferences, "Background top");
    auto* history = find_widget<Fl_Input>(preferences, "Command history");
    if (!first || !history) {
      return fail("preference controls were not found");
    }
    first->value("12 34 56");
    history->value("250");
    if (!preferences.apply() || !applied || changed.background_first.red != 12 ||
        changed.background_first.green != 34 ||
        changed.background_first.blue != 56 || changed.command_history != 250) {
      return fail("valid preferences were not applied: " +
                  preferences.last_error());
    }
    auto* font = find_widget<Fl_Input>(preferences, "Status font (auto or size)");
    auto* terminal_color = find_widget<Fl_Input>(preferences, "Terminal background");
    if (!font || !terminal_color || changed.faceplate_font_size != 0 ||
        changed.terminal_background.red != 255 || changed.terminal_foreground.red != 0) {
      return fail("adaptive status font or terminal color defaults are wrong");
    }
    font->value("18");
    terminal_color->value("30 40 50");
    if (!preferences.apply() || changed.faceplate_font_size != 18 ||
        changed.terminal_background.blue != 50) return fail("appearance preferences failed");
    font->value("auto");
    if (!preferences.apply() || changed.faceplate_font_size != 0)
      return fail("automatic status font was not restored");
    font->value("-5");
    if (preferences.apply()) return fail("negative status font was accepted");
    font->value("auto");
    first->value("300 0 0");
    if (preferences.apply() || preferences.last_error().empty()) {
      return fail("preferences accepted an invalid display color");
    }
  }

  {
    Fl_Double_Window window(700, 420, "View control test");
    window.begin();
    auto* view = new Fl_Cad_View(0, 0, 440, 420, backend);
    auto* controls = new Fl_Ged_View_Controls(440, 0, 260, 420);
    window.end();
    window.show();
    Fl::check();
    const Fl_Cad_View::Color first{12, 34, 56};
    const Fl_Cad_View::Color second{65, 43, 21};
    view->background(first, second);
    if (view->background_first().red != first.red ||
        view->background_second().blue != second.blue) {
      return fail("view display colors were not retained");
    }
    if (!view->initialize() || !view->framebuffer()) {
      return fail("view did not create an embedded framebuffer: " +
                  view->last_error());
    }
    view->resize(0, 0, 420, 400);
    auto resized = view->capture();
    if (!resized || fb_getwidth(view->framebuffer()) != resized->width ||
        fb_getheight(view->framebuffer()) != resized->height) {
      return fail("embedded framebuffer did not track a display resize");
    }
    controls->view(view);
    if (!controls->apply()) {
      return fail("view controls rejected their current values: " +
                  controls->last_error());
    }
    const double initial_scale = view->view()->gv_scale;
    Fl_Cad_View::InputEvent mouse;
    mouse.type = FL_PUSH;
    mouse.button = FL_LEFT_MOUSE;
    mouse.x = mouse.y = 100;
    view->process_input(mouse);
    if (!EQUAL(view->view()->gv_scale, initial_scale * 2.0))
      return fail("MGED left click did not zoom out by two");
    mouse.button = FL_RIGHT_MOUSE;
    view->process_input(mouse);
    if (!EQUAL(view->view()->gv_scale, initial_scale))
      return fail("MGED right click did not zoom in by two");
    mouse.type = FL_MOUSEWHEEL;
    mouse.wheel_y = -1;
    view->process_input(mouse);
    if (!EQUAL(view->view()->gv_scale, initial_scale / 1.1))
      return fail("MGED wheel zoom direction or factor is incorrect");
    for (int button = FL_LEFT_MOUSE; button <= FL_RIGHT_MOUSE; ++button) {
      for (int state : {FL_SHIFT, FL_CTRL, FL_SHIFT | FL_CTRL}) {
        mouse = {};
        mouse.type = FL_PUSH;
        mouse.button = button;
        mouse.buttons = FL_BUTTON(button);
        mouse.state = state;
        mouse.x = mouse.y = 100;
        const double scale = view->view()->gv_scale;
        mat_t rotation, center;
        MAT_COPY(rotation, view->view()->gv_rotation);
        MAT_COPY(center, view->view()->gv_center);
        view->process_input(mouse);
        if (!EQUAL(view->view()->gv_scale, scale))
          return fail("modified mouse press applied an unmodified click zoom");
        mouse.type = FL_DRAG;
        mouse.x += 30;
        mouse.y += 10;
        view->process_input(mouse);
        bool rotated = false, translated = false;
        for (int element = 0; element < 16; ++element) {
          rotated = rotated || !NEAR_EQUAL(rotation[element], view->view()->gv_rotation[element], kViewTolerance);
          translated = translated || !NEAR_EQUAL(center[element], view->view()->gv_center[element], kViewTolerance);
        }
        const bool scaled = !EQUAL(scale, view->view()->gv_scale);
        if (rotated != (state == FL_CTRL) || translated != (state == FL_SHIFT) ||
            scaled != (state == (FL_SHIFT | FL_CTRL)))
          return fail("MGED modifier drag selected the wrong view operation: button=" +
                      std::to_string(button) + " state=" + std::to_string(state) +
                      " rotated=" + std::to_string(rotated) +
                      " translated=" + std::to_string(translated) +
                      " scaled=" + std::to_string(scaled));
        mouse.type = FL_RELEASE;
        view->process_input(mouse);
      }
    }
    Fl_Ged_Interaction interaction(0, 0, 240, 120);
    interaction.ged(gedp);
    interaction.mode(Fl_Ged_Interaction::Mode::QueryRay);
    std::string command;
    interaction.command_handler([&command](const std::string& value) {
      command = value;
      return BRLCAD_OK;
    });
    Fl_Cad_View::InputEvent event;
    event.type = FL_PUSH;
    event.button = FL_MIDDLE_MOUSE;
    event.x = view->w() / 2;
    event.y = view->h() / 2;
    if (!interaction.process_input(*view, event) ||
        command != "vquery_ray 0 0") {
      return fail("view interaction did not issue a center query ray");
    }
    struct bview* original_view = gedp->ged_gvp;
    auto original_callback = gedp->ged_create_vlist_scene_obj_callback;
    gedp->ged_gvp = view->view();
    bv_set_add_view(&gedp->ged_views, view->view());
    gedp->ged_create_vlist_scene_obj_callback = [](void*, struct bv_scene_obj* object) { object->s_flag = UP; };
    view->ged(gedp);
    const char* create_preview[] = {"in", "fltkged_preview.s", "sph", "0", "0", "0", "1"};
    const char* draw_preview[] = {"draw", "fltkged_preview.s"};
    if (ged_exec(gedp, 7, create_preview) != BRLCAD_OK || ged_exec(gedp, 2, draw_preview) != BRLCAD_OK)
      return fail("unable to draw viewport edit fixture");
    struct bv_scene_obj* preview_object = nullptr;
    struct bu_ptbl* displayed = bv_view_objs(view->view(), BV_DB_OBJS);
    for (size_t index = 0; displayed && index < BU_PTBL_LEN(displayed); ++index) {
      auto* object = reinterpret_cast<struct bv_scene_obj*>(BU_PTBL_GET(displayed, index));
      if (!object->s_u_data) continue;
      auto* data = static_cast<struct ged_bv_data*>(object->s_u_data);
      if (data->s_fullpath.fp_len && std::string(DB_FULL_PATH_CUR_DIR(&data->s_fullpath)->d_namep) == "fltkged_preview.s")
        preview_object = object;
    }
    if (!preview_object) return fail("viewport fixture has no displayed scene object");
    const auto minimum_x = [preview_object]() {
      point_t minimum, maximum;
      VSETALL(minimum, INFINITY);
      VSETALL(maximum, -INFINITY);
      size_t length = 0;
      int display_mode = 0;
      bv_vlist_bbox(&preview_object->s_vlist, &minimum, &maximum, &length, &display_mode);
      return minimum[X];
    };
    const double original_minimum = minimum_x();
    double rendered_minimum = original_minimum;
    view->draw_handler([&](struct bview* active) {
      rendered_minimum = minimum_x();
      dm_draw_objs(active, nullptr, nullptr);
    });
    const char* begin_preview[] = {"sed", "fltkged_preview.s"};
    const char* move_preview[] = {"translate", "5", "0", "0"};
    if (ged_edit_session_exec(gedp, 2, begin_preview) != BRLCAD_OK ||
        ged_edit_session_exec(gedp, 4, move_preview) != BRLCAD_OK)
      return fail("unable to move viewport working edit");
    view->invalidate();
    if (!view->capture() || !NEAR_EQUAL(rendered_minimum, original_minimum + 5 * gedp->dbip->dbi_local2base, kViewTolerance))
      return fail("main viewport did not render shared working geometry");
    if (!NEAR_EQUAL(minimum_x(), original_minimum, kViewTolerance))
      return fail("preview changed the persistent display vectors");
    ged_edit_session_end(gedp);
    view->invalidate();
    if (!view->capture() || !NEAR_EQUAL(rendered_minimum, original_minimum, kViewTolerance))
      return fail("reject did not restore the original viewport geometry");
    view->draw_handler({});
    const char* zap_preview[] = {"zap"};
    ged_exec(gedp, 1, zap_preview);
    view->ged(nullptr);
    bv_set_rm_view(&gedp->ged_views, view->view());
    gedp->ged_gvp = original_view;
    gedp->ged_create_vlist_scene_obj_callback = original_callback;
    window.hide();
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
