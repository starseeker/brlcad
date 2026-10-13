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

#include "fltkcad/Fl_Ged_Specialist_Tools.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <FL/Fl_Tabs.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "ged.h"
#include "raytrace.h"
}

namespace {

constexpr int kWindowWidth = 690;
constexpr int kWindowHeight = 420;
constexpr int kMargin = 12;
constexpr int kLabelWidth = 145;
constexpr int kInputHeight = 28;
constexpr int kRowStep = 36;
constexpr int kButtonHeight = 30;
constexpr int kButtonGap = 8;
constexpr std::size_t kMaxFieldBytes = 4096;
constexpr std::size_t kMaxCommandBytes = 4 * kMaxFieldBytes;
constexpr std::uintmax_t kMaxDensityFileBytes = 64U * 1024U * 1024U;
constexpr const char* kMaxDensityFileLabel = "64 MiB";
constexpr double kMaximumNmgAngleDegrees = 180.0;

std::string quote(const std::string& value)
{
  std::string result;
  result.reserve(value.size() + 2);
  result.push_back('"');
  for (const char character : value) {
    if (character == '\\' || character == '"') {
      result.push_back('\\');
    }
    result.push_back(character);
  }
  result.push_back('"');
  return result;
}

bool parse_nonnegative(const char* text, double& value)
{
  if (!text || !*text) {
    return false;
  }
  errno = 0;
  char* end = nullptr;
  value = std::strtod(text, &end);
  return errno != ERANGE && end != text && *end == '\0' &&
         std::isfinite(value) && value >= 0.0;
}

bool parse_rgb(const std::string& text)
{
  if (text == "." || text == "del") {
    return true;
  }
  std::istringstream input(text);
  int values[3]{};
  if (!(input >> values[0] >> values[1] >> values[2])) {
    return false;
  }
  input >> std::ws;
  return input.eof() &&
         std::all_of(std::begin(values), std::end(values),
                     [](int value) { return value >= 0 && value <= 255; });
}

bool valid_density_record(const std::string& record)
{
  const std::size_t first = record.find(',');
  const std::size_t second =
    first == std::string::npos ? first : record.find(',', first + 1);
  if (first == std::string::npos || second == std::string::npos ||
      first == 0 || second == first + 1 || second + 1 == record.size()) {
    return false;
  }

  const std::string id_text = record.substr(0, first);
  const std::string density_text = record.substr(first + 1, second - first - 1);
  char* id_end = nullptr;
  char* density_end = nullptr;
  errno = 0;
  const long id = std::strtol(id_text.c_str(), &id_end, 10);
  const double density = std::strtod(density_text.c_str(), &density_end);
  return errno != ERANGE && id_end && *id_end == '\0' && id >= 0 &&
         density_end && *density_end == '\0' && std::isfinite(density) &&
         density > 0.0;
}

} // namespace

class Fl_Ged_Specialist_Tools::Impl
{
public:
  explicit Impl(Fl_Ged_Specialist_Tools& owner) : owner_(owner)
  {
    owner_.begin();
    tabs_ = new Fl_Tabs(kMargin, kMargin, owner_.w() - 2 * kMargin,
                        owner_.h() - 2 * kMargin);
    tabs_->begin();
    make_bot_tab();
    make_nmg_tab();
    make_material_tab();
    tabs_->end();
    owner_.end();
    owner_.resizable(tabs_);
  }

  void selection(const std::string& name)
  {
    bot_object_->value(name.c_str());
    nmg_object_->value(name.c_str());
    material_object_->value(name.c_str());
    if (!name.empty()) {
      bot_output_->value((name + ".repaired").c_str());
      nmg_output_->value((name + ".simplified").c_str());
    }
  }

  bool inspect_bot()
  {
    const auto object = input(bot_object_);
    return require_object(object, ID_BOT, "BOT") &&
           run("bot check solid " + quote(object), false);
  }

  bool repair_bot()
  {
    const auto source = input(bot_object_);
    const auto output = input(bot_output_);
    return require_object(source, ID_BOT, "BOT") && require_new_name(output) &&
           run("bot repair -o " + quote(output) + " " + quote(source), true);
  }

  bool condense_bot()
  {
    return run_bot_copy("bot_condense");
  }

  bool fuse_bot_vertices()
  {
    return run_bot_copy("bot_vertex_fuse");
  }

  bool fuse_bot_faces()
  {
    return run_bot_copy("bot_face_fuse");
  }

  bool decimate_bot()
  {
    const auto source = input(bot_object_);
    const auto output = input(bot_output_);
    if (!require_object(source, ID_BOT, "BOT") || !require_new_name(output)) {
      return false;
    }

    std::vector<std::string> options;
    const std::pair<const char*, Fl_Input*> fields[] = {
      {"-f", bot_feature_}, {"-c", bot_chord_},
      {"-n", bot_normal_}, {"-e", bot_edge_}};
    for (const auto& field : fields) {
      const std::string value = input(field.second);
      if (value.empty()) {
        continue;
      }
      double parsed = 0.0;
      if (!parse_nonnegative(value.c_str(), parsed)) {
        return fail("BOT decimation tolerances must be finite, nonnegative numbers");
      }
      options.emplace_back(field.first);
      options.push_back(value);
    }
    if (options.empty()) {
      return fail("Specify at least one BOT decimation tolerance");
    }
    if (!input(bot_feature_).empty() && options.size() != 2) {
      return fail("Feature size cannot be combined with the other decimation tolerances");
    }

    std::string command = "bot_decimate";
    for (const auto& option : options) {
      command += " " + quote(option);
    }
    command += " " + quote(output) + " " + quote(source);
    return run(command, true);
  }

  bool synchronize_bot()
  {
    const auto object = input(bot_object_);
    return require_object(object, ID_BOT, "BOT") &&
           run("bot sync " + quote(object), true);
  }

  bool flip_bot()
  {
    const auto object = input(bot_object_);
    return require_object(object, ID_BOT, "BOT") &&
           run("bot flip " + quote(object), true);
  }

  bool fix_nmg_normals()
  {
    const auto object = input(nmg_object_);
    return require_object(object, ID_NMG, "NMG") &&
           run("nmg_fix_normals " + quote(object), true);
  }

  bool collapse_nmg()
  {
    const auto object = input(nmg_object_);
    const auto output = input(nmg_output_);
    double maximum_error = 0.0;
    if (!require_object(object, ID_NMG, "NMG") || !require_new_name(output)) {
      return false;
    }
    if (!parse_nonnegative(nmg_error_->value(), maximum_error)) {
      return fail("NMG collapse error must be a finite, nonnegative number");
    }
    std::string command = "nmg_collapse " + quote(object) + " " +
                          quote(output) + " " + quote(input(nmg_error_));
    const std::string angle = input(nmg_angle_);
    if (!angle.empty()) {
      double parsed = 0.0;
      if (!parse_nonnegative(angle.c_str(), parsed) ||
          parsed > kMaximumNmgAngleDegrees) {
        return fail("NMG minimum angle must be between 0 and 180 degrees");
      }
      command += " " + quote(angle);
    }
    return run(command, true);
  }

  bool simplify_nmg()
  {
    const auto object = input(nmg_object_);
    const auto output = input(nmg_output_);
    if (!require_object(object, ID_NMG, "NMG") || !require_new_name(output)) {
      return false;
    }
    const char* type = nmg_type_->text();
    if (!type || !*type) {
      return fail("Select an NMG simplification target type");
    }
    return run("nmg_simplify " + std::string(type) + " " + quote(output) +
                 " " + quote(object),
               true);
  }

  bool apply_material()
  {
    const auto object = input(material_object_);
    if (!require_object(object, ID_COMBINATION, "combination")) {
      return false;
    }
    const std::string shader = input(material_shader_);
    const std::string color = input(material_color_);
    const std::string inherit = input(material_inherit_);
    if (shader.empty() || shader.size() > kMaxFieldBytes) {
      return fail("Enter a bounded shader value, '.' to retain it, or 'del' to clear it");
    }
    if (!parse_rgb(color)) {
      return fail("Material color must be three values from 0 to 255, '.', or 'del'");
    }
    if (inherit != "y" && inherit != "n" && inherit != ".") {
      return fail("Material inheritance must be y, n, or .");
    }
    std::string color_arguments;
    if (color == "del") {
      color_arguments = quote(color);
    } else if (color == ".") {
      color_arguments = ". . .";
    } else {
      std::istringstream channels(color);
      int red = 0;
      int green = 0;
      int blue = 0;
      channels >> red >> green >> blue;
      color_arguments = std::to_string(red) + " " + std::to_string(green) +
                        " " + std::to_string(blue);
    }
    return run("mater -s " + quote(object) + " " + quote(shader) + " " +
                 color_arguments + " " + quote(inherit), true);
  }

  bool list_densities()
  {
    return require_database() && run("mater -d get", false);
  }

  bool set_density()
  {
    const std::string record = input(density_record_);
    if (!require_database()) {
      return false;
    }
    if (record.empty() || record.size() > kMaxFieldBytes ||
        !valid_density_record(record)) {
      return fail("Density must be entered as a nonnegative id, positive density, and name");
    }
    return run("mater -d set " + quote(record), true);
  }

  bool import_densities(const std::string& path)
  {
    if (!require_database() || !require_density_input(path)) {
      return false;
    }
    if (!run("mater -d validate " + quote(path), false)) {
      return false;
    }
    return run("mater -d import -v " + quote(path), true);
  }

  bool export_densities(const std::string& path)
  {
    return require_database() && require_density_output(path) &&
           run("mater -d export " + quote(path), false);
  }

  bool validate_densities(const std::string& path)
  {
    return require_database() && require_density_input(path) &&
           run("mater -d validate " + quote(path), false);
  }

  bool clear_densities()
  {
    return require_database() && run("mater -d clear", true);
  }

  struct ged* ged_ = nullptr;
  CommandHandler command_handler_;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
  std::string last_error_;

private:
  static std::string input(const Fl_Input* field)
  {
    return field && field->value() ? field->value() : "";
  }

  Fl_Input* add_input(Fl_Group& group, int y, const char* label,
                      const char* initial = "")
  {
    group.begin();
    auto* field = new Fl_Input(group.x() + kMargin + kLabelWidth, y,
                               group.w() - kLabelWidth - 2 * kMargin,
                               kInputHeight, label);
    field->value(initial);
    group.end();
    return field;
  }

  Fl_Button* add_button(Fl_Group& group, int x, int y, int width,
                        const char* label, Fl_Callback callback)
  {
    group.begin();
    auto* button = new Fl_Button(x, y, width, kButtonHeight, label);
    button->callback(callback, this);
    group.end();
    return button;
  }

  void make_bot_tab()
  {
    bot_group_ = new Fl_Group(tabs_->x(), tabs_->y() + kInputHeight,
                              tabs_->w(), tabs_->h() - kInputHeight, "BOT");
    const int left = bot_group_->x() + kMargin;
    int y = bot_group_->y() + kMargin;
    bot_object_ = add_input(*bot_group_, y, "BOT object");
    y += kRowStep;
    bot_output_ = add_input(*bot_group_, y, "Output BOT");
    y += kRowStep;
    bot_feature_ = add_input(*bot_group_, y, "Feature size");
    y += kRowStep;
    bot_chord_ = add_input(*bot_group_, y, "Chord error");
    y += kRowStep;
    bot_normal_ = add_input(*bot_group_, y, "Normal error");
    y += kRowStep;
    bot_edge_ = add_input(*bot_group_, y, "Minimum edge");
    y += kRowStep + 2;
    constexpr int width = 118;
    add_button(*bot_group_, left, y, width, "Inspect", inspect_callback);
    add_button(*bot_group_, left + width + kButtonGap, y, width, "Decimate",
               decimate_callback);
    add_button(*bot_group_, left + 2 * (width + kButtonGap), y, width,
               "Condense", condense_callback);
    add_button(*bot_group_, left + 3 * (width + kButtonGap), y, width,
               "Fuse Vertices", fuse_vertices_callback);
    y += kButtonHeight + kButtonGap;
    add_button(*bot_group_, left, y, width, "Fuse Faces", fuse_faces_callback);
    add_button(*bot_group_, left + width + kButtonGap, y, width, "Sync Normals",
               sync_callback);
    add_button(*bot_group_, left + 2 * (width + kButtonGap), y, width,
               "Flip Normals", flip_callback);
    add_button(*bot_group_, left + 3 * (width + kButtonGap), y, width,
               "Repair Copy", repair_callback);
    bot_group_->end();
  }

  void make_nmg_tab()
  {
    nmg_group_ = new Fl_Group(tabs_->x(), tabs_->y() + kInputHeight,
                              tabs_->w(), tabs_->h() - kInputHeight, "NMG");
    const int left = nmg_group_->x() + kMargin;
    int y = nmg_group_->y() + kMargin;
    nmg_object_ = add_input(*nmg_group_, y, "NMG object");
    y += kRowStep;
    nmg_output_ = add_input(*nmg_group_, y, "Output object");
    y += kRowStep;
    nmg_error_ = add_input(*nmg_group_, y, "Maximum error", "0.01");
    y += kRowStep;
    nmg_angle_ = add_input(*nmg_group_, y, "Minimum angle");
    y += kRowStep;
    nmg_group_->begin();
    nmg_type_ = new Fl_Choice(nmg_group_->x() + kMargin + kLabelWidth, y, 180,
                              kInputHeight, "Simplify as");
    nmg_type_->add("arb");
    nmg_type_->add("tgc");
    nmg_type_->add("poly");
    nmg_type_->value(0);
    nmg_group_->end();
    y += kRowStep + 2;
    constexpr int width = 135;
    add_button(*nmg_group_, left, y, width, "Fix Normals", fix_normals_callback);
    add_button(*nmg_group_, left + width + kButtonGap, y, width, "Collapse Copy",
               collapse_callback);
    add_button(*nmg_group_, left + 2 * (width + kButtonGap), y, width,
               "Simplify Copy", simplify_callback);
    add_button(*nmg_group_, left + 3 * (width + kButtonGap), y, width,
               "Label Faces", label_faces_callback);
    y += kButtonHeight + kButtonGap;
    add_button(*nmg_group_, left, y, width, "Label Vertices",
               label_vertices_callback);
    nmg_group_->end();
  }

  void make_material_tab()
  {
    material_group_ = new Fl_Group(tabs_->x(), tabs_->y() + kInputHeight,
                                   tabs_->w(), tabs_->h() - kInputHeight,
                                   "Materials");
    const int left = material_group_->x() + kMargin;
    int y = material_group_->y() + kMargin;
    material_object_ = add_input(*material_group_, y, "Combination");
    y += kRowStep;
    material_shader_ = add_input(*material_group_, y, "Shader", ".");
    y += kRowStep;
    material_color_ = add_input(*material_group_, y, "RGB", ".");
    y += kRowStep;
    material_inherit_ = add_input(*material_group_, y, "Inherit", ".");
    y += kRowStep;
    density_record_ = add_input(*material_group_, y, "Density record");
    y += kRowStep + 2;
    constexpr int width = 145;
    add_button(*material_group_, left, y, width, "Apply Material",
               apply_material_callback);
    add_button(*material_group_, left + width + kButtonGap, y, width,
               "List Densities", list_densities_callback);
    add_button(*material_group_, left + 2 * (width + kButtonGap), y, width,
               "Set Density", set_density_callback);
    add_button(*material_group_, left + 3 * (width + kButtonGap), y, width,
               "Density Source", density_source_callback);
    y += kButtonHeight + kButtonGap;
    add_button(*material_group_, left, y, width, "Import Densities...",
               import_densities_callback);
    add_button(*material_group_, left + width + kButtonGap, y, width,
               "Export Densities...", export_densities_callback);
    add_button(*material_group_, left + 2 * (width + kButtonGap), y, width,
               "Validate File...", validate_densities_callback);
    add_button(*material_group_, left + 3 * (width + kButtonGap), y, width,
               "Clear Densities", clear_densities_callback);
    material_group_->end();
  }

  bool require_database()
  {
    return ged_ && ged_->dbip ? true : fail("Open a database before using specialist tools");
  }

  bool require_object(const std::string& name, int expected_type,
                      const char* expected_label)
  {
    if (!require_database()) {
      return false;
    }
    if (name.empty() || name.size() > kMaxFieldBytes ||
        name.find('\0') != std::string::npos) {
      return fail(std::string("Enter a bounded ") + expected_label + " object name");
    }
    const struct directory* object =
      db_lookup(ged_->dbip, name.c_str(), LOOKUP_QUIET);
    if (!object || object->d_minor_type != expected_type) {
      return fail(name + " is not a " + expected_label + " object");
    }
    return true;
  }

  bool require_new_name(const std::string& name)
  {
    if (!require_database()) {
      return false;
    }
    if (name.empty() || name.size() > kMaxFieldBytes ||
        name.find('\0') != std::string::npos) {
      return fail("Enter a bounded output object name");
    }
    if (db_lookup(ged_->dbip, name.c_str(), LOOKUP_QUIET)) {
      return fail("The output object already exists: " + name);
    }
    return true;
  }

  bool require_density_input(const std::string& path)
  {
    if (!valid_path(path)) {
      return fail("Select a bounded density file path");
    }
    std::error_code error;
    const std::filesystem::path file(path);
    if (!std::filesystem::is_regular_file(file, error) || error) {
      return fail("The selected density file is missing or is not a regular file");
    }
    const std::uintmax_t size = std::filesystem::file_size(file, error);
    if (error) {
      return fail("Unable to determine the selected density file size");
    }
    if (size > kMaxDensityFileBytes) {
      return fail(std::string("The selected density file exceeds the ") +
                  kMaxDensityFileLabel + " safety limit");
    }
    return true;
  }

  bool require_density_output(const std::string& path)
  {
    if (!valid_path(path)) {
      return fail("Select a bounded density output path");
    }
    std::error_code error;
    const std::filesystem::path file(path);
    if (std::filesystem::exists(file, error) || error) {
      return fail(error ? "Unable to inspect the density output path"
                        : "The density output file already exists");
    }
    const std::filesystem::path parent = file.parent_path();
    if (!parent.empty() &&
        (!std::filesystem::is_directory(parent, error) || error)) {
      return fail("The density output directory does not exist");
    }
    return true;
  }

  static bool valid_path(const std::string& path)
  {
    return !path.empty() && path.size() <= kMaxFieldBytes &&
           path.find('\0') == std::string::npos;
  }

  bool run_bot_copy(const char* command)
  {
    const auto source = input(bot_object_);
    const auto output = input(bot_output_);
    return require_object(source, ID_BOT, "BOT") && require_new_name(output) &&
           run(std::string(command) + " " + quote(output) + " " +
                 quote(source),
               true);
  }

  bool run(const std::string& command, bool changed)
  {
    if (!command_handler_) {
      return fail("No command handler is connected to the specialist tools");
    }
    if (command.size() > kMaxCommandBytes) {
      return fail("The generated specialist command exceeds its size limit");
    }
    if (command_handler_(command) != BRLCAD_OK) {
      return fail("BRL-CAD rejected the specialist operation");
    }
    last_error_.clear();
    if (changed && changed_handler_) {
      changed_handler_();
    }
    return true;
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return false;
  }

  template <bool (Impl::*Operation)()>
  static void operation_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!(self.*Operation)()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void inspect_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::inspect_bot>(w, d);
  }

  static void repair_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::repair_bot>(w, d);
  }

  static void condense_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::condense_bot>(w, d);
  }

  static void fuse_vertices_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::fuse_bot_vertices>(w, d);
  }

  static void fuse_faces_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::fuse_bot_faces>(w, d);
  }

  static void decimate_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::decimate_bot>(w, d);
  }

  static void sync_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::synchronize_bot>(w, d);
  }

  static void flip_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::flip_bot>(w, d);
  }

  static void fix_normals_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::fix_nmg_normals>(w, d);
  }

  static void collapse_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::collapse_nmg>(w, d);
  }

  static void simplify_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::simplify_nmg>(w, d);
  }

  static void apply_material_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::apply_material>(w, d);
  }

  static void list_densities_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::list_densities>(w, d);
  }

  static void set_density_callback(Fl_Widget* w, void* d)
  {
    operation_callback<&Impl::set_density>(w, d);
  }

  enum class DensityFileAction
  {
    Import,
    Export,
    Validate
  };

  void choose_density_file(DensityFileAction action)
  {
    Fl_Native_File_Chooser chooser;
    chooser.filter("Density file\t*.density");
    if (action == DensityFileAction::Export) {
      chooser.title("Export density table");
      chooser.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
      chooser.options(Fl_Native_File_Chooser::SAVEAS_CONFIRM |
                      Fl_Native_File_Chooser::USE_FILTER_EXT);
    } else {
      chooser.title(action == DensityFileAction::Import
                      ? "Import density table"
                      : "Validate density file");
      chooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
    }
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }

    bool succeeded = false;
    const std::string path = chooser.filename();
    if (action == DensityFileAction::Import) {
      if (fl_choice("Import replaces the current database density table. Continue?",
                    "Cancel", "Import", nullptr) != 1) {
        return;
      }
      succeeded = import_densities(path);
    } else if (action == DensityFileAction::Export) {
      succeeded = export_densities(path);
    } else {
      succeeded = validate_densities(path);
    }
    if (!succeeded) {
      fl_alert("%s", last_error_.c_str());
    }
  }

  static void import_densities_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->choose_density_file(DensityFileAction::Import);
  }

  static void export_densities_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->choose_density_file(DensityFileAction::Export);
  }

  static void validate_densities_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->choose_density_file(DensityFileAction::Validate);
  }

  static void clear_densities_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (fl_choice("Remove all density records from the database?", "Cancel",
                  "Clear", nullptr) != 1) {
      return;
    }
    if (!self.clear_densities()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void density_source_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if ((!self.require_database() || !self.run("mater -d source", false))) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void label_faces_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const auto object = input(self.nmg_object_);
    if ((!self.require_object(object, ID_NMG, "NMG") ||
         !self.run("labelface " + quote(object), false))) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void label_vertices_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const auto object = input(self.nmg_object_);
    if ((!self.require_object(object, ID_NMG, "NMG") ||
         !self.run("labelvert " + quote(object), false))) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  Fl_Ged_Specialist_Tools& owner_;
  Fl_Tabs* tabs_ = nullptr;
  Fl_Group* bot_group_ = nullptr;
  Fl_Input* bot_object_ = nullptr;
  Fl_Input* bot_output_ = nullptr;
  Fl_Input* bot_feature_ = nullptr;
  Fl_Input* bot_chord_ = nullptr;
  Fl_Input* bot_normal_ = nullptr;
  Fl_Input* bot_edge_ = nullptr;
  Fl_Group* nmg_group_ = nullptr;
  Fl_Input* nmg_object_ = nullptr;
  Fl_Input* nmg_output_ = nullptr;
  Fl_Input* nmg_error_ = nullptr;
  Fl_Input* nmg_angle_ = nullptr;
  Fl_Choice* nmg_type_ = nullptr;
  Fl_Group* material_group_ = nullptr;
  Fl_Input* material_object_ = nullptr;
  Fl_Input* material_shader_ = nullptr;
  Fl_Input* material_color_ = nullptr;
  Fl_Input* material_inherit_ = nullptr;
  Fl_Input* density_record_ = nullptr;
};

Fl_Ged_Specialist_Tools::Fl_Ged_Specialist_Tools()
  : Fl_Double_Window(kWindowWidth, kWindowHeight, "Specialist Geometry Tools"),
    impl_(std::make_unique<Impl>(*this))
{
}

Fl_Ged_Specialist_Tools::~Fl_Ged_Specialist_Tools() = default;
void Fl_Ged_Specialist_Tools::ged(struct ged* value) { impl_->ged_ = value; }
void Fl_Ged_Specialist_Tools::selection(const std::string& object_name)
{
  impl_->selection(object_name);
}
void Fl_Ged_Specialist_Tools::command_handler(CommandHandler handler)
{
  impl_->command_handler_ = std::move(handler);
}
void Fl_Ged_Specialist_Tools::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler_ = std::move(handler);
}
void Fl_Ged_Specialist_Tools::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}
bool Fl_Ged_Specialist_Tools::inspect_bot() { return impl_->inspect_bot(); }
bool Fl_Ged_Specialist_Tools::repair_bot() { return impl_->repair_bot(); }
bool Fl_Ged_Specialist_Tools::condense_bot() { return impl_->condense_bot(); }
bool Fl_Ged_Specialist_Tools::fuse_bot_vertices()
{
  return impl_->fuse_bot_vertices();
}
bool Fl_Ged_Specialist_Tools::fuse_bot_faces()
{
  return impl_->fuse_bot_faces();
}
bool Fl_Ged_Specialist_Tools::decimate_bot() { return impl_->decimate_bot(); }
bool Fl_Ged_Specialist_Tools::synchronize_bot()
{
  return impl_->synchronize_bot();
}
bool Fl_Ged_Specialist_Tools::flip_bot() { return impl_->flip_bot(); }
bool Fl_Ged_Specialist_Tools::fix_nmg_normals()
{
  return impl_->fix_nmg_normals();
}
bool Fl_Ged_Specialist_Tools::collapse_nmg() { return impl_->collapse_nmg(); }
bool Fl_Ged_Specialist_Tools::simplify_nmg() { return impl_->simplify_nmg(); }
bool Fl_Ged_Specialist_Tools::apply_material()
{
  return impl_->apply_material();
}
bool Fl_Ged_Specialist_Tools::list_densities()
{
  return impl_->list_densities();
}
bool Fl_Ged_Specialist_Tools::set_density() { return impl_->set_density(); }
bool Fl_Ged_Specialist_Tools::import_densities(const std::string& path)
{
  return impl_->import_densities(path);
}
bool Fl_Ged_Specialist_Tools::export_densities(const std::string& path)
{
  return impl_->export_densities(path);
}
bool Fl_Ged_Specialist_Tools::validate_densities(const std::string& path)
{
  return impl_->validate_densities(path);
}
bool Fl_Ged_Specialist_Tools::clear_densities()
{
  return impl_->clear_densities();
}
const std::string& Fl_Ged_Specialist_Tools::last_error() const noexcept
{
  return impl_->last_error_;
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
