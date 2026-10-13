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

#include "fltkcad/Fl_Ged_Primitive_Editor.H"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Scroll.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

extern "C" {
#include "bg.h"
#include "bu/vls.h"
#include "bv/vlist.h"
#include "dm.h"
#include "ged.h"
#include "ged/edit.h"
#include "raytrace.h"
#include "rt/edit.h"
#include "rt/functab.h"
}

namespace {

constexpr int kInitialWidth = 1000;
constexpr int kInitialHeight = 650;
constexpr int kMargin = 8;
constexpr int kControlsWidth = 350;
constexpr int kRowHeight = 28;
constexpr int kLabelWidth = 132;
constexpr int kButtonHeight = 30;
constexpr int kMaxCommands = 512;
constexpr int kMaxParameters = 64;
constexpr std::size_t kMaxTextBytes = RT_EDIT_MAXSTR_LEN - 1;
constexpr std::size_t kMaxLogBytes = 64U * 1024U;
constexpr int kPreviewRed = 255;
constexpr int kPreviewGreen = 215;
constexpr int kPreviewBlue = 0;

bool list_parameter(const struct rt_edit_param_desc& parameter)
{
  return parameter.type == RT_EDIT_PARAM_SCALAR_LIST ||
    parameter.type == RT_EDIT_PARAM_INTEGER_LIST;
}

int parameter_width(const struct rt_edit_param_desc& parameter)
{
  switch (parameter.type) {
    case RT_EDIT_PARAM_POINT2:
    case RT_EDIT_PARAM_VECTOR2: return 2;
    case RT_EDIT_PARAM_POINT:
    case RT_EDIT_PARAM_VECTOR:
    case RT_EDIT_PARAM_COLOR: return 3;
    case RT_EDIT_PARAM_MATRIX: return 16;
    case RT_EDIT_PARAM_SCALAR_LIST:
    case RT_EDIT_PARAM_INTEGER_LIST: return parameter.nenum;
    case RT_EDIT_PARAM_STRING: return 0;
    default: return 1;
  }
}

bool parse_values(const char* text, int minimum, int maximum,
                  std::vector<fastf_t>& values)
{
  values.clear();
  const char* cursor = text;
  while (*cursor) {
    while (*cursor && std::isspace(static_cast<unsigned char>(*cursor))) {
      ++cursor;
    }
    if (!*cursor) {
      break;
    }
    if (values.size() == static_cast<std::size_t>(maximum)) {
      return false;
    }
    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(cursor, &end);
    if (cursor == end || errno == ERANGE || !std::isfinite(value)) {
      return false;
    }
    values.push_back(value);
    cursor = end;
  }
  return values.size() >= static_cast<std::size_t>(minimum);
}

std::string format_values(const fastf_t* values, int count)
{
  std::ostringstream output;
  output << std::setprecision(12);
  for (int index = 0; index < count; ++index) {
    if (index != 0) {
      output << ' ';
    }
    output << values[index];
  }
  return output.str();
}

} // namespace

class Fl_Ged_Primitive_Editor::Impl
{
public:
  struct Field
  {
    const struct rt_edit_param_desc* description = nullptr;
    Fl_Box* label = nullptr;
    Fl_Input* input = nullptr;
    Fl_Check_Button* check = nullptr;
    Fl_Choice* choice = nullptr;
  };

  Impl(Fl_Ged_Primitive_Editor& owner, Fl_Cad_View::Backend backend)
    : owner_(owner)
  {
    owner_.begin();
    heading_ = new Fl_Box(kMargin, kMargin, kControlsWidth - 2 * kMargin,
                          kRowHeight, "No primitive loaded");
    heading_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    command_ = new Fl_Choice(kMargin + kLabelWidth, kMargin + kRowHeight,
                             kControlsWidth - kLabelWidth - 2 * kMargin,
                             kRowHeight, "Operation");
    command_->align(FL_ALIGN_LEFT);
    command_->callback(command_changed, this);

    fields_scroll_ = new Fl_Scroll(kMargin, kMargin + 2 * kRowHeight,
                                   kControlsWidth - 2 * kMargin,
                                   owner_.h() - 3 * kRowHeight -
                                     kButtonHeight - 3 * kMargin);
    fields_scroll_->begin();
    fields_group_ = new Fl_Group(fields_scroll_->x(), fields_scroll_->y(),
                                 fields_scroll_->w(), fields_scroll_->h());
    fields_group_->end();
    fields_scroll_->end();

    apply_ = new Fl_Button(kMargin, owner_.h() - kButtonHeight - kMargin,
                           76, kButtonHeight, "Apply");
    revert_ = new Fl_Button(kMargin + 82,
                            owner_.h() - kButtonHeight - kMargin, 76,
                            kButtonHeight, "Undo");
    save_ = new Fl_Button(kMargin + 164,
                          owner_.h() - kButtonHeight - kMargin, 76,
                          kButtonHeight, "Save");
    close_ = new Fl_Button(kMargin + 246,
                           owner_.h() - kButtonHeight - kMargin, 88,
                           kButtonHeight, "Close");
    apply_->callback(apply_callback, this);
    revert_->callback(revert_callback, this);
    save_->callback(save_callback, this);
    close_->callback(close_callback, this);

    preview_ = new Fl_Cad_View(kControlsWidth, 0,
                               owner_.w() - kControlsWidth, owner_.h(),
                               backend);
    preview_->draw_handler([this](struct bview* view) { draw(view); });
    preview_->error_handler(
      [this](const std::string& message) { report_error(message); });
    owner_.end();
    owner_.resizable(preview_);
    owner_.callback(close_callback, this);
    set_enabled(false);
  }

  ~Impl()
  {
    clear_fields();
    if (active()) ged_edit_interactive_end(ged_);
  }

  bool open(struct ged* gedp, const std::string& object_name)
  {
    if (dirty()) {
      return fail("Save or discard the current primitive changes before opening another object");
    }
    reset();
    if (!gedp || !gedp->dbip) {
      return fail("Open a database before editing a primitive");
    }
    if (object_name.empty() || object_name.find('\0') != std::string::npos) {
      return fail("Select a valid primitive name");
    }

    struct db_full_path path = DB_FULL_PATH_INIT_ZERO;
    if (db_string_to_path(&path, gedp->dbip, object_name.c_str()) < 0 || !path.fp_len) {
      db_free_full_path(&path);
      return fail("The selected primitive no longer exists");
    }
    struct directory* candidate = DB_FULL_PATH_CUR_DIR(&path);
    if (candidate->d_flags & RT_DIR_COMB) {
      db_free_full_path(&path);
      return fail("Use the combination editor for combination objects");
    }
    if (ged_edit_interactive_begin(gedp, nullptr, &path, GED_EDIT_PRIMITIVE,
                              path.fp_len - 1, &tolerance_, gedp->ged_gvp) != BRLCAD_OK) {
      db_free_full_path(&path);
      return fail(bu_vls_cstr(gedp->ged_result_str));
    }
    char* selected = db_path_to_string(&path);
    std::string session_path = selected;
    bu_free(selected, "panel edit selection");
    db_free_full_path(&path);

    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    const struct rt_edit_prim_desc* description = nullptr;
    if (ged_edit_session_default_get(gedp, session_path.c_str(),
            GED_EDIT_TARGET_OBJECT, 0, &session) != GED_EDIT_OK ||
        ged_edit_session_descriptor_get(gedp, session, &description) != GED_EDIT_OK ||
        description->ncmd <= 0 || description->ncmd > kMaxCommands) {
      ged_edit_interactive_end(gedp);
      return fail("This primitive has no bounded structured edit interface");
    }

    ged_ = gedp;
    edit_ = ged_edit_interactive_edit(gedp);
    description_ = description;
    object_name_ = object_name;
    session_path_ = std::move(session_path);
    (void)ged_edit_session_checkpoint(ged_, session);

    command_->clear();
    for (int index = 0; index < description_->ncmd; ++index) {
      const char* label = description_->cmds[index].label;
      command_->add(label ? label : "Unnamed operation");
    }
    command_->value(0);
    heading_->copy_label((std::string(description_->prim_label
                                       ? description_->prim_label
                                       : "Primitive") +
                          ": " + object_name_)
                           .c_str());
    owner_.copy_label(("Primitive Editor - " + object_name_).c_str());
    set_enabled(true);
    configure_command();
    fit_preview();
    preview_->invalidate();
    last_error_.clear();
    return true;
  }

  bool apply()
  {
    if (!active() || !description_) {
      return fail("No primitive is open for editing");
    }
    const int selected = command_->value();
    if (selected < 0 || selected >= description_->ncmd) {
      return fail("Select a valid edit operation");
    }
    const struct rt_edit_cmd_desc& command = description_->cmds[selected];
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    if (!session_get(session)) {
      return fail("The edit session is stale");
    }

    fastf_t values[RT_EDIT_MAXPARA] = {};
    const char* strings[RT_EDIT_MAXSTR] = {};
    struct ged_edit_command_input input = {command.cmd_id, values, 0, strings, 0};
    std::vector<fastf_t> parsed;
    for (const auto& field : fields_) {
      const auto& parameter = *field.description;
      if (parameter.type == RT_EDIT_PARAM_STRING) {
        const char* value = field.input->value();
        if (std::strlen(value) > kMaxTextBytes) {
          return fail(std::string(parameter.label) + " is too long");
        }
        if (parameter.index < 0 || parameter.index >= RT_EDIT_MAXSTR) {
          return fail("Primitive string parameter index is invalid");
        }
        strings[parameter.index] = value;
        input.string_count = std::max(input.string_count,
          static_cast<std::size_t>(parameter.index) + 1);
        continue;
      }

      const int width = parameter_width(parameter);
      if (parameter.index < 0 || width <= 0 ||
          parameter.index > RT_EDIT_MAXPARA - width) {
        return fail("Primitive numeric parameter index is invalid");
      }
      parsed.clear();
      if (parameter.type == RT_EDIT_PARAM_BOOLEAN) {
        parsed.push_back(field.check->value() ? 1.0 : 0.0);
      } else if (parameter.type == RT_EDIT_PARAM_ENUM) {
        const int choice = field.choice->value();
        if (choice < 0 || choice >= parameter.nenum || !parameter.enum_ids) {
          return fail(std::string("Select a value for ") + parameter.label);
        }
        parsed.push_back(parameter.enum_ids[choice]);
      } else {
        const int maximum = list_parameter(parameter) ? RT_EDIT_MAXPARA - parameter.index : width;
        if (!parse_values(field.input->value(), width, maximum, parsed)) {
          return fail(std::string("Invalid numeric values for ") + parameter.label);
        }
      }

      std::copy(parsed.begin(), parsed.end(), values + parameter.index);
      input.value_count = std::max(input.value_count, parameter.index + parsed.size());
    }

    if (ged_edit_session_checkpoint(ged_, session) != GED_EDIT_OK) {
      return fail("Unable to create an undo checkpoint");
    }
    bu_vls_trunc(ged_->ged_result_str, 0);
    const enum ged_edit_status status = ged_edit_session_apply(ged_, session, &input);
    if (status != GED_EDIT_OK) {
      if (status == GED_EDIT_REJECTED) {
        const auto message = edit_log("BRL-CAD rejected the requested edit");
        (void)ged_edit_session_revert(ged_, session);
        return fail(message);
      }
      return fail("The requested edit has invalid parameter values");
    }
    if (changed_handler_) changed_handler_();
    configure_command();
    preview_->invalidate();
    update_buttons();
    last_error_.clear();
    return true;
  }

  bool revert()
  {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    if (!session_get(session) || ged_edit_session_revert(ged_, session) != GED_EDIT_OK) {
      return fail("There is no edit checkpoint to restore");
    }
    if (changed_handler_) changed_handler_();
    configure_command();
    preview_->invalidate();
    update_buttons();
    last_error_.clear();
    return true;
  }

  bool save()
  {
    if (!active()) return fail("No primitive is open for saving");
    if (ged_edit_interactive_apply(ged_, 1) != BRLCAD_OK)
      return fail(bu_vls_cstr(ged_->ged_result_str));
    edit_ = ged_edit_interactive_edit(ged_);
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    if (session_get(session)) {
      (void)ged_edit_session_checkpoint(ged_, session);
    }
    configure_command();
    update_buttons();
    preview_->invalidate();
    if (changed_handler_) changed_handler_();
    last_error_.clear();
    return true;
  }

  void resize(int width, int height)
  {
    const int controls_width = std::min(kControlsWidth, std::max(1, width / 2));
    heading_->resize(kMargin, kMargin, controls_width - 2 * kMargin,
                     kRowHeight);
    command_->resize(kMargin + kLabelWidth, kMargin + kRowHeight,
                     controls_width - kLabelWidth - 2 * kMargin, kRowHeight);
    fields_scroll_->resize(kMargin, kMargin + 2 * kRowHeight,
                           controls_width - 2 * kMargin,
                           std::max(1, height - 3 * kRowHeight -
                                         kButtonHeight - 3 * kMargin));
    apply_->position(kMargin, height - kButtonHeight - kMargin);
    revert_->position(kMargin + 82, height - kButtonHeight - kMargin);
    save_->position(kMargin + 164, height - kButtonHeight - kMargin);
    close_->position(kMargin + 246, height - kButtonHeight - kMargin);
    preview_->resize(controls_width, 0, std::max(1, width - controls_width),
                     height);
    layout_fields();
  }

  void changed_handler(ChangedHandler handler)
  {
    changed_handler_ = std::move(handler);
  }

  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

  bool dirty() const noexcept
  {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    struct ged_edit_session_info info = {};
    if (!active()) {
      return false;
    }
    /* A stale session can still hold edits that the user must discard. */
    if (!session_get(session) ||
        ged_edit_session_info_get(ged_, session, &info) != GED_EDIT_OK) {
      return true;
    }
    return info.dirty;
  }
  const std::string& object_name() const noexcept { return object_name_; }
  const std::string& last_error() const noexcept { return last_error_; }

private:
  bool active() const noexcept
  {
    if (!ged_ || !edit_ || ged_edit_interactive_edit(ged_) != edit_ ||
        ged_edit_interactive_mode(ged_) != GED_EDIT_PRIMITIVE) return false;
    return ged_edit_interactive_matches(ged_, session_path_.c_str(), GED_EDIT_PRIMITIVE);
  }

  bool session_get(ged_edit_session_ref& session) const noexcept
  {
    return active() && ged_edit_session_default_get(ged_, session_path_.c_str(),
      GED_EDIT_TARGET_OBJECT, 0, &session) == GED_EDIT_OK;
  }

  void reset()
  {
    clear_fields();
    command_->clear();
    if (active()) ged_edit_interactive_end(ged_);
    edit_ = nullptr;
    session_path_.clear();
    ged_ = nullptr;
    description_ = nullptr;
    object_name_.clear();
    heading_->copy_label("No primitive loaded");
    set_enabled(false);
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* widgets[] = {command_, apply_, revert_, save_};
    for (auto* widget : widgets) {
      if (enabled) {
        widget->activate();
      } else {
        widget->deactivate();
      }
    }
    update_buttons();
  }

  void update_buttons()
  {
    if (dirty()) {
      save_->activate();
    } else {
      save_->deactivate();
    }
    if (edit_) {
      revert_->activate();
    } else {
      revert_->deactivate();
    }
  }

  void clear_fields()
  {
    for (auto& field : fields_) {
      delete field.label;
      delete field.input;
      delete field.check;
      delete field.choice;
    }
    fields_.clear();
  }

  void configure_command()
  {
    clear_fields();
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    if (!description_ || !session_get(session)) {
      return;
    }
    const int selected = command_->value();
    if (selected < 0 || selected >= description_->ncmd) {
      return;
    }
    const auto& command = description_->cmds[selected];
    if (command.nparam < 0 || command.nparam > kMaxParameters) {
      report_error("The edit operation exceeds the field limit");
      return;
    }

    struct rt_edit_cmd_values current;
    if (ged_edit_session_command_values_get(ged_, session, command.cmd_id, &current) == GED_EDIT_ERROR) {
      report_error("Could not read the current edit values");
    }
    fields_group_->begin();
    for (int index = 0; index < command.nparam; ++index) {
      const auto& parameter = command.params[index];
      Field field;
      field.description = &parameter;
      field.label = new Fl_Box(0, 0, kLabelWidth, kRowHeight,
                               parameter.label ? parameter.label : "Value");
      field.label->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
      if (parameter.type == RT_EDIT_PARAM_BOOLEAN) {
        field.check = new Fl_Check_Button(0, 0, 40, kRowHeight);
        if (parameter.index >= 0 && parameter.index < RT_EDIT_MAXPARA) {
          field.check->value(!NEAR_ZERO(current.values[parameter.index], SMALL_FASTF));
        }
      } else if (parameter.type == RT_EDIT_PARAM_ENUM) {
        field.choice = new Fl_Choice(0, 0, 100, kRowHeight);
        int selected_enum = 0;
        for (int option = 0; option < parameter.nenum; ++option) {
          field.choice->add(parameter.enum_labels && parameter.enum_labels[option]
                              ? parameter.enum_labels[option]
                              : "Unnamed");
          if (parameter.enum_ids && parameter.index >= 0 &&
              parameter.index < RT_EDIT_MAXPARA &&
              parameter.enum_ids[option] ==
                static_cast<int>(current.values[parameter.index])) {
            selected_enum = option;
          }
        }
        field.choice->value(selected_enum);
      } else {
        field.input = new Fl_Input(0, 0, 100, kRowHeight);
        if (parameter.type == RT_EDIT_PARAM_STRING) {
          if (parameter.index >= 0 && parameter.index < RT_EDIT_MAXSTR &&
              current.string_valid[parameter.index]) {
            field.input->value(current.strings[parameter.index]);
          }
        } else {
          int width = parameter_width(parameter);
          if (list_parameter(parameter) && parameter.index >= 0 &&
              current.value_count > static_cast<std::size_t>(parameter.index)) {
            width = static_cast<int>(current.value_count) - parameter.index;
          }
          if (parameter.index >= 0 && width > 0 &&
              parameter.index + width <= RT_EDIT_MAXPARA) {
            field.input->value(
              format_values(current.values + parameter.index, width).c_str());
          }
        }
      }
      fields_.push_back(field);
    }
    fields_group_->end();
    layout_fields();
    fields_scroll_->redraw();
  }

  void layout_fields()
  {
    const int x = fields_scroll_->x();
    const int width = std::max(1, fields_scroll_->w() - 20);
    const int value_x = x + kLabelWidth;
    const int value_width = std::max(1, width - kLabelWidth);
    int y = fields_scroll_->y();
    for (auto& field : fields_) {
      field.label->resize(x, y, kLabelWidth, kRowHeight);
      if (field.input) {
        field.input->resize(value_x, y, value_width, kRowHeight);
      }
      if (field.check) {
        field.check->resize(value_x, y, value_width, kRowHeight);
      }
      if (field.choice) {
        field.choice->resize(value_x, y, value_width, kRowHeight);
      }
      y += kRowHeight + 4;
    }
    fields_group_->resize(x, fields_scroll_->y(), width,
                          std::max(fields_scroll_->h(), y - fields_scroll_->y()));
  }

  void fit_preview()
  {
    if (!active() || !preview_->view()) {
      return;
    }
    const int type = edit_->es_int.idb_type;
    point_t minimum;
    point_t maximum;
    if (!OBJ[type].ft_bbox ||
        OBJ[type].ft_bbox(&edit_->es_int, &minimum, &maximum, &tolerance_) < 0) {
      return;
    }
    point_t center;
    VADD2SCALE(center, minimum, maximum, 0.5);
    MAT_DELTAS_VEC_NEG(preview_->view()->gv_center, center);
    const double diagonal = DIST_PNT_PNT(minimum, maximum);
    preview_->view()->gv_scale = std::max(1.0, diagonal * 0.65);
    preview_->view()->gv_size = 2.0 * preview_->view()->gv_scale;
    preview_->view()->gv_isize = 1.0 / preview_->view()->gv_size;
    bv_update(preview_->view());
  }

  void draw(struct bview* view)
  {
    if (!active() || !view || !view->dmp) {
      return;
    }
    struct bu_list vlist;
    BU_LIST_INIT(&vlist);
    struct bg_tess_tol tessellation;
    tessellation.magic = BG_TESS_TOL_MAGIC;
    tessellation.abs = 0.0;
    tessellation.rel = 0.01;
    tessellation.norm = 0.0;
    const int type = edit_->es_int.idb_type;
    if (OBJ[type].ft_plot &&
        OBJ[type].ft_plot(&vlist, &edit_->es_int, &tessellation, &tolerance_,
                          view) == 0) {
      auto* display = static_cast<struct dm*>(view->dmp);
      dm_set_fg(display, kPreviewRed, kPreviewGreen, kPreviewBlue, 1, 1.0);
      dm_draw_vlist(display, reinterpret_cast<struct bv_vlist*>(&vlist));
    }
    bv_vlist_cleanup(&vlist);
  }

  std::string edit_log(const char* fallback) const
  {
    if (ged_ && bu_vls_strlen(ged_->ged_result_str) > 0) {
      const std::size_t length = bu_vls_strlen(ged_->ged_result_str);
      const std::size_t accepted = std::min(length, kMaxLogBytes);
      std::string result(bu_vls_cstr(ged_->ged_result_str), accepted);
      if (accepted != length) {
        result += "\n[edit diagnostic truncated]";
      }
      return result;
    }
    return fallback;
  }

  bool fail(std::string message)
  {
    report_error(std::move(message));
    return false;
  }

  void report_error(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
  }

  bool request_close()
  {
    if (!dirty()) {
      reset();
      if (changed_handler_) changed_handler_();
      owner_.hide();
      return true;
    }
    const int decision = fl_choice(
      "Save changes to %s before closing?", "Cancel", "Discard", "Save",
      object_name_.c_str());
    if (decision == 0) {
      return false;
    }
    if (decision == 2 && !save()) {
      fl_alert("%s", last_error_.c_str());
      return false;
    }
    reset();
    if (changed_handler_) changed_handler_();
    owner_.hide();
    return true;
  }

  static void command_changed(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->configure_command();
  }

  static void apply_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.apply()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void revert_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.revert()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void save_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.save()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void close_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->request_close();
  }

  Fl_Ged_Primitive_Editor& owner_;
  Fl_Box* heading_ = nullptr;
  Fl_Choice* command_ = nullptr;
  Fl_Scroll* fields_scroll_ = nullptr;
  Fl_Group* fields_group_ = nullptr;
  Fl_Button* apply_ = nullptr;
  Fl_Button* revert_ = nullptr;
  Fl_Button* save_ = nullptr;
  Fl_Button* close_ = nullptr;
  Fl_Cad_View* preview_ = nullptr;
  std::vector<Field> fields_;
  struct ged* ged_ = nullptr;
  struct rt_edit* edit_ = nullptr;
  const struct rt_edit_prim_desc* description_ = nullptr;
  struct bn_tol tolerance_ = BN_TOL_INIT_TOL;
  std::string session_path_;
  std::string object_name_;
  std::string last_error_;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_Primitive_Editor::Fl_Ged_Primitive_Editor(
  Fl_Cad_View::Backend backend)
  : Fl_Double_Window(kInitialWidth, kInitialHeight, "Primitive Editor"),
    impl_(std::make_unique<Impl>(*this, backend))
{
}

Fl_Ged_Primitive_Editor::~Fl_Ged_Primitive_Editor() = default;

bool Fl_Ged_Primitive_Editor::open(struct ged* gedp,
                                   const std::string& object_name)
{
  return impl_->open(gedp, object_name);
}

bool Fl_Ged_Primitive_Editor::apply() { return impl_->apply(); }
bool Fl_Ged_Primitive_Editor::revert() { return impl_->revert(); }
bool Fl_Ged_Primitive_Editor::save() { return impl_->save(); }
bool Fl_Ged_Primitive_Editor::dirty() const noexcept { return impl_->dirty(); }
const std::string& Fl_Ged_Primitive_Editor::object_name() const noexcept
{
  return impl_->object_name();
}
const std::string& Fl_Ged_Primitive_Editor::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Primitive_Editor::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler(std::move(handler));
}
void Fl_Ged_Primitive_Editor::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
}
void Fl_Ged_Primitive_Editor::resize(int x, int y, int width, int height)
{
  Fl_Double_Window::resize(x, y, width, height);
  if (impl_) {
    impl_->resize(width, height);
  }
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
