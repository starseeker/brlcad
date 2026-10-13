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

#include "fltkcad/Fl_Ged_Interaction.H"

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bu/vls.h"
#include "bv/util.h"
#include "ged.h"
#include "ged/edit.h"
}

namespace {

constexpr int kMargin = 8;
constexpr int kRowHeight = 27;
constexpr int kRowGap = 5;
constexpr int kLabelWidth = 54;
constexpr int kPickRadiusPixels = 8;
constexpr double kMinimumZoomFraction = 0.01;
constexpr double kDragRotationDegreesPerPixel = 0.35;
constexpr double kDragScalePerPixel = 0.01;

std::string number(double value)
{
  std::ostringstream output;
  output << std::setprecision(15) << value;
  return output.str();
}

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

std::string object_from_path(const std::string& path)
{
  const std::size_t separator = path.find_last_of('/');
  std::string object =
    separator == std::string::npos ? path : path.substr(separator + 1);
  const std::size_t instance = object.find('@');
  if (instance != std::string::npos) {
    object.erase(instance);
  }
  return object;
}

} // namespace

class Fl_Ged_Interaction::Impl
{
public:
  Impl(Fl_Ged_Interaction& owner, int x, int y, int width, int)
    : owner_(owner)
  {
    owner_.begin();
    mode_ = new Fl_Choice(0, 0, 1, kRowHeight, "Mode");
    mode_->add("Navigate|Select|Query ray|Rectangle select|Rectangle zoom|"
               "Move object|Rotate object|Scale object");
    mode_->value(0);
    mode_->callback(mode_callback, this);
    selection_ = new Fl_Box(0, 0, 1, kRowHeight, "No object selected");
    selection_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE | FL_ALIGN_CLIP);
    cancel_ = new Fl_Button(0, 0, 1, kRowHeight, "Cancel gesture");
    cancel_->callback(cancel_callback, this);
    owner_.end();
    layout(x, y, width);
  }

  void ged(struct ged* value)
  {
    ged_ = value;
    if (!ged_) {
      clear_gesture();
    }
  }

  void selection(std::string path, std::string object)
  {
    selected_path_ = std::move(path);
    selected_object_ = std::move(object);
    selection_->copy_label(selected_path_.empty()
                             ? "No object selected"
                             : ("Selected: " + selected_path_).c_str());
  }

  Mode mode() const noexcept
  {
    return static_cast<Mode>(std::clamp(mode_->value(), 0, 7));
  }

  void mode(Mode value)
  {
    clear_gesture();
    mode_->value(static_cast<int>(value));
  }

  bool handle(Fl_Cad_View& widget, const Fl_Cad_View::InputEvent& event)
  {
    if ((event.state & (FL_SHIFT | FL_CTRL | FL_ALT)) != 0) return false;
    if (mode() == Mode::Navigate) {
      if (event.type == FL_PUSH && event.button == FL_MIDDLE_MOUSE &&
          ged_ && ged_->dbip && widget.view()) {
        return pick(widget, event.x, event.y);
      }
      return false;
    }
    if (!ged_ || !ged_->dbip || !widget.view()) {
      return fail("Open a database before using an interaction tool");
    }
    if (event.button != FL_MIDDLE_MOUSE && event.type == FL_PUSH) {
      return false;
    }

    switch (event.type) {
      case FL_PUSH:
        begin(widget, event.x, event.y);
        if (mode() == Mode::Select) {
          return pick(widget, event.x, event.y);
        }
        if (mode() == Mode::QueryRay) {
          return query(widget, event.x, event.y);
        }
        return true;
      case FL_DRAG:
        if (!dragging_) {
          return false;
        }
        update(widget, event.x, event.y);
        return true;
      case FL_RELEASE:
        if (!dragging_) {
          return false;
        }
        update(widget, event.x, event.y);
        return finish(widget);
      case FL_KEYDOWN:
        if (event.key == FL_Escape) {
          clear_gesture();
          widget.invalidate();
          return true;
        }
        return false;
      default:
        return false;
    }
  }

  const std::string& last_error() const noexcept { return last_error_; }

  void resize(int x, int y, int width, int height)
  {
    static_cast<void>(height);
    layout(x, y, width);
  }

  Fl_Ged_Interaction::CommandHandler command_handler_;
  Fl_Ged_Interaction::SelectionHandler selection_handler_;
  Fl_Ged_Interaction::ChangedHandler changed_handler_;
  Fl_Ged_Interaction::ErrorHandler error_handler_;

private:
  static std::pair<double, double> normalized(const Fl_Cad_View& widget,
                                               int x, int y)
  {
    const double width = std::max(1, widget.w());
    const double height = std::max(1, widget.h());
    return {2.0 * static_cast<double>(x) / width - 1.0,
            1.0 - 2.0 * static_cast<double>(y) / height};
  }

  void begin(Fl_Cad_View& widget, int x, int y)
  {
    dragging_ = true;
    start_x_ = end_x_ = x;
    start_y_ = end_y_ = y;
    if (uses_rectangle()) {
      update_rectangle(widget);
    }
  }

  void update(Fl_Cad_View& widget, int x, int y)
  {
    end_x_ = std::clamp(x, 0, std::max(0, widget.w() - 1));
    end_y_ = std::clamp(y, 0, std::max(0, widget.h() - 1));
    if (uses_rectangle()) {
      update_rectangle(widget);
      widget.invalidate();
    }
  }

  bool finish(Fl_Cad_View& widget)
  {
    const Mode current = mode();
    bool result = true;
    if (current == Mode::RectangleSelect) {
      result = pick_rectangle(widget);
    } else if (current == Mode::RectangleZoom) {
      result = zoom_rectangle(widget);
    } else if (current == Mode::TranslateObject ||
               current == Mode::RotateObject ||
               current == Mode::ScaleObject) {
      result = transform(widget, current);
    }
    clear_gesture();
    widget.invalidate();
    return result;
  }

  bool pick(Fl_Cad_View& widget, int x, int y)
  {
    const auto point = normalized(widget, x, y);
    const double radius =
      2.0 * kPickRadiusPixels / std::max(1, std::min(widget.w(), widget.h()));
    const std::string command = "select -p " + number(point.first) + " " +
                                number(point.second) + " " + number(radius);
    const bool selected = select_from_command(command);
    clear_gesture();
    widget.invalidate();
    return selected;
  }

  bool pick_rectangle(Fl_Cad_View& widget)
  {
    const auto start = normalized(widget, start_x_, start_y_);
    const auto end = normalized(widget, end_x_, end_y_);
    const std::string command =
      "select -p " + number(start.first) + " " + number(start.second) + " " +
      number(end.first - start.first) + " " + number(end.second - start.second);
    return select_from_command(command);
  }

  bool select_from_command(const std::string& command)
  {
    if (!ged_) {
      return fail("GED is not initialized");
    }
    std::vector<std::string> arguments;
    std::istringstream input(command);
    for (std::string argument; input >> argument;) {
      arguments.push_back(std::move(argument));
    }
    std::vector<const char*> argv;
    argv.reserve(arguments.size());
    for (const auto& argument : arguments) {
      argv.push_back(argument.c_str());
    }
    bu_vls_trunc(ged_->ged_result_str, 0);
    ged_->ged_gvp = active_view_;
    const int status = ged_exec(ged_, static_cast<int>(argv.size()), argv.data());
    const std::string output = bu_vls_cstr(ged_->ged_result_str);
    bu_vls_trunc(ged_->ged_result_str, 0);
    if (status != BRLCAD_OK) {
      return fail(output.empty() ? "Geometry selection failed" : output);
    }

    std::istringstream paths(output);
    std::string path;
    if (!std::getline(paths, path) || path.empty()) {
      return fail("No displayed geometry was found at that location");
    }
    const std::string object = object_from_path(path);
    if (object.empty() || !db_lookup(ged_->dbip, object.c_str(), LOOKUP_QUIET)) {
      return fail("The selected display path does not resolve to a database object");
    }
    if (!selected_object_.empty()) {
      const char* clear[] = {"illum", "-n", selected_object_.c_str()};
      bu_vls_trunc(ged_->ged_result_str, 0);
      ged_exec(ged_, 3, clear);
    }
    const char* illuminate[] = {"illum", object.c_str()};
    bu_vls_trunc(ged_->ged_result_str, 0);
    ged_exec(ged_, 2, illuminate);
    bu_vls_trunc(ged_->ged_result_str, 0);
    selection(path, object);
    if (selection_handler_) {
      selection_handler_({path, object});
    }
    if (changed_handler_) {
      changed_handler_();
    }
    last_error_.clear();
    return true;
  }

  bool query(Fl_Cad_View& widget, int x, int y)
  {
    const auto point = normalized(widget, x, y);
    const long view_x = std::lround(point.first * BV_MAX);
    const long view_y = std::lround(point.second * BV_MAX);
    clear_gesture();
    if (!command_handler_) {
      return fail("No command handler is configured for query rays");
    }
    command_handler_("vquery_ray " + std::to_string(view_x) + " " +
                     std::to_string(view_y));
    return true;
  }

  bool zoom_rectangle(Fl_Cad_View& widget)
  {
    struct bview* view = widget.view();
    const auto start = normalized(widget, start_x_, start_y_);
    const auto end = normalized(widget, end_x_, end_y_);
    const double width = std::abs(end.first - start.first);
    const double height = std::abs(end.second - start.second);
    const double fraction = std::max(width, height) * 0.5;
    if (fraction < kMinimumZoomFraction) {
      return fail("The zoom rectangle is too small");
    }
    point_t view_center;
    VSET(view_center, (start.first + end.first) * 0.5,
         (start.second + end.second) * 0.5, 0.0);
    point_t model_center;
    MAT4X3PNT(model_center, view->gv_view2model, view_center);
    bv_view_set_center_vec(view, model_center);
    bv_view_set_size(view, view->gv_size * fraction);
    bv_update(view);
    if (changed_handler_) {
      changed_handler_();
    }
    last_error_.clear();
    return true;
  }

  bool transform(Fl_Cad_View& widget, Mode operation)
  {
    if (selected_path_.empty()) {
      return fail("Select an object before starting a transform");
    }
    if (!command_handler_) {
      return fail("No command handler is configured for transforms");
    }
    const int dx = end_x_ - start_x_;
    const int dy = end_y_ - start_y_;
    if (dx == 0 && dy == 0) {
      return true;
    }

    ged_->ged_gvp = widget.view();
    if (ged_edit_session_mode(ged_) == GED_EDIT_VIEW) {
      if (command_handler_("oed " + quote(selected_path_)) != BRLCAD_OK) return true;
    } else {
      if (!ged_edit_session_matches(ged_, selected_path_.c_str(), GED_EDIT_INSTANCE))
        return fail("Accept or reject the active edit before moving another selection");
    }
    std::string command;
    if (operation == Mode::TranslateObject) {
      const auto start = normalized(widget, start_x_, start_y_);
      const auto end = normalized(widget, end_x_, end_y_);
      point_t start_view;
      point_t end_view;
      VSET(start_view, start.first, start.second, 0.0);
      VSET(end_view, end.first, end.second, 0.0);
      point_t start_model;
      point_t end_model;
      MAT4X3PNT(start_model, widget.view()->gv_view2model, start_view);
      MAT4X3PNT(end_model, widget.view()->gv_view2model, end_view);
      vect_t delta;
      VSUB2(delta, end_model, start_model);
      const double base2local = widget.view()->gv_base2local > 0.0
                                  ? widget.view()->gv_base2local
                                  : 1.0;
      VSCALE(delta, delta, base2local);
      command = "tra -m " + number(delta[X]) +
                " " + number(delta[Y]) + " " + number(delta[Z]);
    } else if (operation == Mode::RotateObject) {
      const double angle = (dx - dy) * kDragRotationDegreesPerPixel;
      command = "rot 0 0 " + number(angle);
    } else {
      const double factor = std::exp((dx - dy) * kDragScalePerPixel);
      command = "sca " + number(factor);
    }
    if (command_handler_(command) == BRLCAD_OK) {
      if (changed_handler_) {
        changed_handler_();
      }
      last_error_.clear();
    }
    return true;
  }

  bool uses_rectangle() const noexcept
  {
    return mode() == Mode::RectangleSelect || mode() == Mode::RectangleZoom;
  }

  void update_rectangle(Fl_Cad_View& widget)
  {
    struct bview_settings* settings = widget.view()->gv_s
                                        ? widget.view()->gv_s
                                        : &widget.view()->gv_ls;
    auto& rectangle = settings->gv_rect;
    const auto start = normalized(widget, start_x_, start_y_);
    const auto end = normalized(widget, end_x_, end_y_);
    rectangle.active = 1;
    rectangle.draw = 1;
    rectangle.pos[0] = start_x_;
    rectangle.pos[1] = start_y_;
    rectangle.dim[0] = end_x_ - start_x_;
    rectangle.dim[1] = end_y_ - start_y_;
    rectangle.x = start.first;
    rectangle.y = start.second;
    rectangle.width = end.first - start.first;
    rectangle.height = end.second - start.second;
    rectangle.cdim[0] = widget.w();
    rectangle.cdim[1] = widget.h();
    rectangle.aspect = static_cast<double>(widget.w()) / std::max(1, widget.h());
  }

  void clear_gesture()
  {
    dragging_ = false;
    if (active_view_) {
      struct bview_settings* settings = active_view_->gv_s
                                          ? active_view_->gv_s
                                          : &active_view_->gv_ls;
      settings->gv_rect.active = 0;
      settings->gv_rect.draw = 0;
    }
    active_view_ = nullptr;
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return true;
  }

  void layout(int x, int y, int width)
  {
    const int content_width = std::max(1, width - 2 * kMargin);
    mode_->resize(x + kMargin + kLabelWidth, y + kMargin,
                  std::max(1, content_width - kLabelWidth), kRowHeight);
    selection_->resize(x + kMargin, y + kMargin + kRowHeight + kRowGap,
                       content_width, kRowHeight);
    cancel_->resize(x + kMargin, y + kMargin + 2 * (kRowHeight + kRowGap),
                    content_width, kRowHeight);
  }

  static void mode_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->clear_gesture();
  }

  static void cancel_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->clear_gesture();
  }

  Fl_Ged_Interaction& owner_;
  struct ged* ged_ = nullptr;
  Fl_Choice* mode_ = nullptr;
  Fl_Box* selection_ = nullptr;
  Fl_Button* cancel_ = nullptr;
  struct bview* active_view_ = nullptr;
  std::string selected_path_;
  std::string selected_object_;
  bool dragging_ = false;
  int start_x_ = 0;
  int start_y_ = 0;
  int end_x_ = 0;
  int end_y_ = 0;
  std::string last_error_;

  void set_active_view(struct bview* view) { active_view_ = view; }

  friend class Fl_Ged_Interaction;
};

Fl_Ged_Interaction::Fl_Ged_Interaction(int x, int y, int width, int height,
                                       const char* label)
  : Fl_Group(x, y, width, height, label),
    impl_(std::make_unique<Impl>(*this, x, y, width, height))
{
}

Fl_Ged_Interaction::~Fl_Ged_Interaction() = default;
void Fl_Ged_Interaction::ged(struct ged* value) { impl_->ged(value); }
void Fl_Ged_Interaction::selection(std::string path, std::string object)
{
  impl_->selection(std::move(path), std::move(object));
}
Fl_Ged_Interaction::Mode Fl_Ged_Interaction::mode() const noexcept
{
  return impl_->mode();
}
void Fl_Ged_Interaction::mode(Mode value) { impl_->mode(value); }
bool Fl_Ged_Interaction::process_input(
  Fl_Cad_View& view, const Fl_Cad_View::InputEvent& event)
{
  impl_->set_active_view(view.view());
  return impl_->handle(view, event);
}
const std::string& Fl_Ged_Interaction::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Interaction::command_handler(CommandHandler handler)
{
  impl_->command_handler_ = std::move(handler);
}
void Fl_Ged_Interaction::selection_handler(SelectionHandler handler)
{
  impl_->selection_handler_ = std::move(handler);
}
void Fl_Ged_Interaction::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler_ = std::move(handler);
}
void Fl_Ged_Interaction::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}
void Fl_Ged_Interaction::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  impl_->resize(x, y, width, height);
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
