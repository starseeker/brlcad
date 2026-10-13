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

#include "fltkcad/Fl_Ged_View_Controls.H"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Scroll.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bv/adc.h"
#include "bv/util.h"
}

namespace {

constexpr int kMargin = 8;
constexpr int kRowHeight = 25;
constexpr int kRowGap = 3;
constexpr int kLabelWidth = 78;
constexpr int kScrollBarWidth = 18;
constexpr int kButtonGap = 6;
constexpr int kButtonWidth = 78;
constexpr int kInputLimit = 256;
constexpr double kMaximumPerspective = 179.0;

bool parse_values(const char* text, double* values, std::size_t count)
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  for (std::size_t index = 0; index < count; ++index) {
    if (!(input >> values[index]) || !std::isfinite(values[index])) {
      return false;
    }
  }
  std::string trailing;
  return !(input >> trailing);
}

std::string format_values(const double* values, std::size_t count)
{
  std::ostringstream output;
  output << std::setprecision(8);
  for (std::size_t index = 0; index < count; ++index) {
    if (index) {
      output << ' ';
    }
    output << values[index];
  }
  return output.str();
}

struct ToggleBinding
{
  const char* label;
  Fl_Check_Button* widget = nullptr;
};

} // namespace

class Fl_Ged_View_Controls::Impl
{
public:
  Impl(Fl_Ged_View_Controls& owner, int x, int y, int width, int height)
    : owner_(owner)
  {
    owner_.begin();
    scroll_ = new Fl_Scroll(x, y, width, height);
    scroll_->begin();
    heading_ = new Fl_Box(0, 0, 1, kRowHeight, "Active View");
    heading_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    aet_ = make_input("AET");
    center_ = make_input("Center");
    size_ = make_input("Size");
    perspective_ = make_input("Perspective");
    grid_anchor_ = make_input("Grid anchor");
    grid_spacing_ = make_input("Grid spacing");
    grid_major_ = make_input("Grid major");
    adc_position_ = make_input("ADC position");
    adc_angles_ = make_input("ADC angles");
    adc_distance_ = make_input("ADC distance");

    const std::array<const char*, 16> labels{
      "Adaptive CSG", "Adaptive mesh", "Angle/distance cursor",
      "Center dot",   "Grid",          "Snap to grid",
      "Model axes",   "View axes",     "Scale",
      "Z clipping",   "Parameters",    "Parameter size",
      "Parameter center", "Parameter azimuth", "Parameter elevation",
      "Parameter twist"};
    toggles_.reserve(labels.size());
    for (const char* label : labels) {
      ToggleBinding binding;
      binding.label = label;
      binding.widget = new Fl_Check_Button(0, 0, 1, kRowHeight, label);
      toggles_.push_back(binding);
    }
    fps_ = new Fl_Check_Button(0, 0, 1, kRowHeight, "Parameter FPS");
    framebuffer_ = new Fl_Choice(0, 0, 1, kRowHeight, "Framebuffer");
    framebuffer_->add("Off|Overlay|Underlay");
    framebuffer_->value(0);
    apply_ = new Fl_Button(0, 0, kButtonWidth, kRowHeight, "Apply");
    refresh_ = new Fl_Button(0, 0, kButtonWidth, kRowHeight, "Refresh");
    apply_->callback(apply_callback, this);
    refresh_->callback(refresh_callback, this);
    scroll_->end();
    owner_.end();
    owner_.resizable(scroll_);
    layout(x, y, width, height);
    set_enabled(false);
  }

  void view(Fl_Cad_View* view)
  {
    view_ = view;
    set_enabled(view_ && view_->view());
    refresh();
  }

  Fl_Cad_View* view() const noexcept { return view_; }

  bool apply()
  {
    struct bview* value = view_ ? view_->view() : nullptr;
    if (!value) {
      return fail("No active view is available");
    }

    double aet[3]{};
    double center[3]{};
    double local_size = 0.0;
    double perspective = 0.0;
    double grid_anchor[3]{};
    double grid_spacing[2]{};
    double grid_major[2]{};
    double adc_position[3]{};
    double adc_angles[2]{};
    double adc_distance = 0.0;
    if (!parse_values(aet_->value(), aet, 3)) {
      return fail("AET must contain three finite numbers");
    }
    if (!parse_values(center_->value(), center, 3)) {
      return fail("Center must contain three finite numbers");
    }
    if (!parse_values(size_->value(), &local_size, 1) || local_size <= 0.0) {
      return fail("View size must be a positive finite number");
    }
    if (!parse_values(perspective_->value(), &perspective, 1) ||
        perspective < 0.0 || perspective > kMaximumPerspective) {
      return fail("Perspective must be between 0 and 179 degrees");
    }
    if (!parse_values(grid_anchor_->value(), grid_anchor, 3) ||
        !parse_values(grid_spacing_->value(), grid_spacing, 2) ||
        grid_spacing[0] <= 0.0 || grid_spacing[1] <= 0.0 ||
        !parse_values(grid_major_->value(), grid_major, 2) ||
        grid_major[0] < 1.0 || grid_major[1] < 1.0 ||
        !EQUAL(std::trunc(grid_major[0]), grid_major[0]) ||
        !EQUAL(std::trunc(grid_major[1]), grid_major[1])) {
      return fail("Grid anchor needs three values; spacing must be positive and major intervals must be positive integers");
    }
    if (!parse_values(adc_position_->value(), adc_position, 3) ||
        !parse_values(adc_angles_->value(), adc_angles, 2) ||
        !parse_values(adc_distance_->value(), &adc_distance, 1) ||
        adc_distance < 0.0) {
      return fail("ADC position, angles, or non-negative distance is invalid");
    }

    const double base2local = valid_scale(value->gv_base2local);
    vect_t model_aet;
    VSET(model_aet, aet[0], aet[1], aet[2]);
    point_t model_center;
    VSET(model_center, center[0] / base2local, center[1] / base2local,
         center[2] / base2local);
    bv_view_set_aet(value, model_aet);
    bv_view_set_center_vec(value, model_center);
    bv_view_set_size(value, local_size / base2local);
    bv_view_set_perspective(value, perspective);

    struct bview_settings* settings =
      value->gv_s ? value->gv_s : &value->gv_ls;
    settings->adaptive_plot_csg = checked(0);
    settings->adaptive_plot_mesh = checked(1);
    settings->gv_adc.draw = checked(2);
    settings->gv_center_dot.gos_draw = checked(3);
    settings->gv_grid.draw = checked(4);
    settings->gv_grid.snap = checked(5);
    settings->gv_model_axes.draw = checked(6);
    settings->gv_view_axes.draw = checked(7);
    settings->gv_view_scale.gos_draw = checked(8);
    settings->gv_zclip = checked(9);
    settings->gv_view_params.draw = checked(10);
    settings->gv_view_params.draw_size = checked(11);
    settings->gv_view_params.draw_center = checked(12);
    settings->gv_view_params.draw_az = checked(13);
    settings->gv_view_params.draw_el = checked(14);
    settings->gv_view_params.draw_tw = checked(15);
    settings->gv_view_params.draw_fps = fps_->value() ? 1 : 0;
    settings->gv_fb_mode = std::clamp(framebuffer_->value(), 0, 2);
    VSET(settings->gv_grid.anchor, grid_anchor[0] / base2local,
         grid_anchor[1] / base2local, grid_anchor[2] / base2local);
    settings->gv_grid.res_h = grid_spacing[0] / base2local;
    settings->gv_grid.res_v = grid_spacing[1] / base2local;
    settings->gv_grid.res_major_h = static_cast<int>(grid_major[0]);
    settings->gv_grid.res_major_v = static_cast<int>(grid_major[1]);
    VSET(settings->gv_adc.pos_model, adc_position[0] / base2local,
         adc_position[1] / base2local, adc_position[2] / base2local);
    settings->gv_adc.a1 = adc_angles[0];
    settings->gv_adc.a2 = adc_angles[1];
    settings->gv_adc.dv_a1 =
      static_cast<int>((1.0 - settings->gv_adc.a1 / 45.0) * BV_MAX);
    settings->gv_adc.dv_a2 =
      static_cast<int>((1.0 - settings->gv_adc.a2 / 45.0) * BV_MAX);
    settings->gv_adc.dst =
      adc_distance / (valid_scale(value->gv_scale) * base2local);
    settings->gv_adc.dv_dist = static_cast<int>(
      (settings->gv_adc.dst / M_SQRT1_2 - 1.0) * BV_MAX);
    adc_model_to_adc_view(&settings->gv_adc, value->gv_model2view, BV_MAX);

    bv_update(value);
    view_->invalidate();
    last_error_.clear();
    if (changed_handler_) {
      changed_handler_();
    }
    return true;
  }

  void refresh()
  {
    struct bview* value = view_ ? view_->view() : nullptr;
    if (!value) {
      aet_->value("");
      center_->value("");
      size_->value("");
      perspective_->value("");
      grid_anchor_->value("");
      grid_spacing_->value("");
      grid_major_->value("");
      adc_position_->value("");
      adc_angles_->value("");
      adc_distance_->value("");
      return;
    }

    vect_t current_aet;
    point_t current_center;
    bv_view_get_aet(value, current_aet);
    bv_view_get_center_vec(value, current_center);
    const double base2local = valid_scale(value->gv_base2local);
    const double aet_values[3]{current_aet[0], current_aet[1], current_aet[2]};
    const double center_values[3]{current_center[0] * base2local,
                                  current_center[1] * base2local,
                                  current_center[2] * base2local};
    const double size_value = bv_view_get_size(value) * base2local;
    const double perspective_value = bv_view_get_perspective(value);
    aet_->value(format_values(aet_values, 3).c_str());
    center_->value(format_values(center_values, 3).c_str());
    size_->value(format_values(&size_value, 1).c_str());
    perspective_->value(format_values(&perspective_value, 1).c_str());

    const struct bview_settings* settings =
      value->gv_s ? value->gv_s : &value->gv_ls;
    set_checked(0, settings->adaptive_plot_csg);
    set_checked(1, settings->adaptive_plot_mesh);
    set_checked(2, settings->gv_adc.draw);
    set_checked(3, settings->gv_center_dot.gos_draw);
    set_checked(4, settings->gv_grid.draw);
    set_checked(5, settings->gv_grid.snap);
    set_checked(6, settings->gv_model_axes.draw);
    set_checked(7, settings->gv_view_axes.draw);
    set_checked(8, settings->gv_view_scale.gos_draw);
    set_checked(9, settings->gv_zclip);
    set_checked(10, settings->gv_view_params.draw);
    set_checked(11, settings->gv_view_params.draw_size);
    set_checked(12, settings->gv_view_params.draw_center);
    set_checked(13, settings->gv_view_params.draw_az);
    set_checked(14, settings->gv_view_params.draw_el);
    set_checked(15, settings->gv_view_params.draw_tw);
    fps_->value(settings->gv_view_params.draw_fps ? 1 : 0);
    framebuffer_->value(std::clamp(settings->gv_fb_mode, 0, 2));
    const double grid_anchor[3]{settings->gv_grid.anchor[0] * base2local,
                                settings->gv_grid.anchor[1] * base2local,
                                settings->gv_grid.anchor[2] * base2local};
    const double grid_spacing[2]{settings->gv_grid.res_h * base2local,
                                 settings->gv_grid.res_v * base2local};
    const double grid_major[2]{
      static_cast<double>(settings->gv_grid.res_major_h),
      static_cast<double>(settings->gv_grid.res_major_v)};
    const double adc_position[3]{settings->gv_adc.pos_model[0] * base2local,
                                 settings->gv_adc.pos_model[1] * base2local,
                                 settings->gv_adc.pos_model[2] * base2local};
    const double adc_angles[2]{settings->gv_adc.a1, settings->gv_adc.a2};
    const double adc_distance =
      settings->gv_adc.dst * valid_scale(value->gv_scale) * base2local;
    grid_anchor_->value(format_values(grid_anchor, 3).c_str());
    grid_spacing_->value(format_values(grid_spacing, 2).c_str());
    grid_major_->value(format_values(grid_major, 2).c_str());
    adc_position_->value(format_values(adc_position, 3).c_str());
    adc_angles_->value(format_values(adc_angles, 2).c_str());
    adc_distance_->value(format_values(&adc_distance, 1).c_str());
    last_error_.clear();
  }

  void resize(int x, int y, int width, int height)
  {
    layout(x, y, width, height);
  }

  const std::string& last_error() const noexcept { return last_error_; }
  void changed_handler(ChangedHandler handler)
  {
    changed_handler_ = std::move(handler);
  }
  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

private:
  Fl_Input* make_input(const char* label)
  {
    auto* input = new Fl_Input(0, 0, 1, kRowHeight, label);
    input->align(FL_ALIGN_LEFT);
    input->maximum_size(kInputLimit);
    return input;
  }

  static double valid_scale(double scale)
  {
    return std::isfinite(scale) && scale > 0.0 ? scale : 1.0;
  }

  int checked(std::size_t index) const
  {
    return toggles_[index].widget->value() ? 1 : 0;
  }

  void set_checked(std::size_t index, int value)
  {
    toggles_[index].widget->value(value ? 1 : 0);
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* fixed[] = {aet_, center_, size_, perspective_, grid_anchor_,
                          grid_spacing_, grid_major_, adc_position_,
                          adc_angles_, adc_distance_, fps_, framebuffer_,
                          apply_, refresh_};
    for (auto* widget : fixed) {
      enabled ? widget->activate() : widget->deactivate();
    }
    for (auto& toggle : toggles_) {
      enabled ? toggle.widget->activate() : toggle.widget->deactivate();
    }
  }

  void layout(int x, int y, int width, int height)
  {
    scroll_->resize(x, y, width, height);
    const int content_width = std::max(1, width - 2 * kMargin - kScrollBarWidth);
    int row_y = y + kMargin;
    heading_->resize(x + kMargin, row_y, content_width, kRowHeight);
    row_y += kRowHeight + kRowGap;
    Fl_Input* inputs[] = {aet_, center_, size_, perspective_, grid_anchor_,
                          grid_spacing_, grid_major_, adc_position_,
                          adc_angles_, adc_distance_};
    for (auto* input : inputs) {
      input->resize(x + kMargin + kLabelWidth, row_y,
                    std::max(1, content_width - kLabelWidth), kRowHeight);
      row_y += kRowHeight + kRowGap;
    }
    row_y += kRowGap;
    for (auto& toggle : toggles_) {
      toggle.widget->resize(x + kMargin, row_y, content_width, kRowHeight);
      row_y += kRowHeight + kRowGap;
    }
    fps_->resize(x + kMargin, row_y, content_width, kRowHeight);
    row_y += kRowHeight + kRowGap;
    framebuffer_->resize(x + kMargin + kLabelWidth, row_y,
                         std::max(1, content_width - kLabelWidth), kRowHeight);
    row_y += kRowHeight + 2 * kRowGap;
    const int buttons_width = 2 * kButtonWidth + kButtonGap;
    const int buttons_x = x + kMargin + std::max(0, content_width - buttons_width);
    apply_->resize(buttons_x, row_y, kButtonWidth, kRowHeight);
    refresh_->resize(buttons_x + kButtonWidth + kButtonGap, row_y,
                     kButtonWidth, kRowHeight);
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return false;
  }

  static void apply_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.apply()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void refresh_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->refresh();
  }

  Fl_Ged_View_Controls& owner_;
  Fl_Scroll* scroll_ = nullptr;
  Fl_Box* heading_ = nullptr;
  Fl_Input* aet_ = nullptr;
  Fl_Input* center_ = nullptr;
  Fl_Input* size_ = nullptr;
  Fl_Input* perspective_ = nullptr;
  Fl_Input* grid_anchor_ = nullptr;
  Fl_Input* grid_spacing_ = nullptr;
  Fl_Input* grid_major_ = nullptr;
  Fl_Input* adc_position_ = nullptr;
  Fl_Input* adc_angles_ = nullptr;
  Fl_Input* adc_distance_ = nullptr;
  std::vector<ToggleBinding> toggles_;
  Fl_Check_Button* fps_ = nullptr;
  Fl_Choice* framebuffer_ = nullptr;
  Fl_Button* apply_ = nullptr;
  Fl_Button* refresh_ = nullptr;
  Fl_Cad_View* view_ = nullptr;
  std::string last_error_;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_View_Controls::Fl_Ged_View_Controls(int x, int y, int width,
                                           int height, const char* label)
  : Fl_Group(x, y, width, height, label),
    impl_(std::make_unique<Impl>(*this, x, y, width, height))
{
}

Fl_Ged_View_Controls::~Fl_Ged_View_Controls() = default;
void Fl_Ged_View_Controls::view(Fl_Cad_View* view) { impl_->view(view); }
Fl_Cad_View* Fl_Ged_View_Controls::view() const noexcept
{
  return impl_->view();
}
bool Fl_Ged_View_Controls::apply() { return impl_->apply(); }
void Fl_Ged_View_Controls::refresh() { impl_->refresh(); }
const std::string& Fl_Ged_View_Controls::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_View_Controls::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler(std::move(handler));
}
void Fl_Ged_View_Controls::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
}
void Fl_Ged_View_Controls::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  if (impl_) {
    impl_->resize(x, y, width, height);
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
