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

#include "fltkcad/Fl_Ged_Transform_Editor.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>

#include <cmath>
#include <sstream>
#include <string>
#include <utility>

extern "C" {
#include "bu/defines.h"
#include "ged/edit.h"
}

namespace {

constexpr int kWindowWidth = 540;
constexpr int kWindowHeight = 270;
constexpr int kMargin = 12;
constexpr int kRowHeight = 28;
constexpr int kGap = 7;
constexpr int kLabelWidth = 112;

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

bool valid_numbers(const char* text, int expected, bool positive = false)
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  double value = 0.0;
  for (int index = 0; index < expected; ++index) {
    if (!(input >> value) || !std::isfinite(value) ||
        (positive && value <= 0.0)) {
      return false;
    }
  }
  std::string trailing;
  return !(input >> trailing);
}

} // namespace

class Fl_Ged_Transform_Editor::Impl
{
public:
  explicit Impl(Fl_Ged_Transform_Editor& owner) : owner_(owner)
  {
    owner_.begin();
    path_ = input("Path", 0);
    target_ = new Fl_Choice(kMargin + kLabelWidth,
                            kMargin + kRowHeight + kGap,
                            kWindowWidth - 2 * kMargin - kLabelWidth,
                            kRowHeight, "Target");
    target_->add("Primitive definition|Selected instance matrix");
    target_->value(0);
    translation_ = input("Translation", 2);
    translation_->value("0 0 0");
    rotation_ = input("Rotation XYZ", 3);
    rotation_->value("0 0 0");
    scale_ = input("Uniform scale", 4);
    scale_->value("1");

    const int button_y = kMargin + 5 * (kRowHeight + kGap);
    translate_ = new Fl_Button(kMargin, button_y, 116, kRowHeight, "Translate");
    rotate_ = new Fl_Button(kMargin + 124, button_y, 116, kRowHeight, "Rotate");
    apply_scale_ = new Fl_Button(kMargin + 248, button_y, 116, kRowHeight,
                                 "Scale");
    reset_ = new Fl_Button(kMargin + 372, button_y, 144, kRowHeight,
                           "Reset Edit");
    for (auto* button : {translate_, rotate_, apply_scale_, reset_}) {
      button->callback(button_callback, this);
    }
    const int finish_y = button_y + kRowHeight + kGap;
    accept_ = new Fl_Button(kMargin, finish_y, 116, kRowHeight, "Accept");
    reject_ = new Fl_Button(kMargin + 124, finish_y, 116, kRowHeight, "Reject");
    accept_->callback(button_callback, this);
    reject_->callback(button_callback, this);
    owner_.end();
    owner_.resizable(path_);
  }

  void selection(std::string path)
  {
    path_->value(path.c_str());
    while (!path.empty() && path.front() == '/') {
      path.erase(path.begin());
    }
    target_->value(path.find('/') == std::string::npos ? 0 : 1);
  }

  void button(Fl_Widget* widget)
  {
    std::string path = path_->value() ? path_->value() : "";
    if (path.empty()) {
      fail("Select an object or enter a database path first");
      return;
    }
    const bool member = target_->value() == 1;
    while (member && !path.empty() && path.front() == '/') {
      path.erase(path.begin());
    }
    if (member && path.find('/') == std::string::npos) {
      fail("A member matrix requires a combination/member path");
      return;
    }

    if (widget == accept_ || widget == reject_ || widget == reset_) {
      run(widget == accept_ ? "accept" : widget == reject_ ? "reject" : "sed_reset");
      return;
    }
    const char* value = widget == translate_ ? translation_->value() :
      widget == rotate_ ? rotation_->value() : scale_->value();
    const bool scaling = widget == apply_scale_;
    if (!valid_numbers(value, scaling ? 1 : 3, scaling)) {
      fail(scaling ? "Scale must be one finite number greater than zero" : "Enter three finite numeric values");
      return;
    }
    if (!ged_) {
      fail("No GED context is configured for editing");
      return;
    }
    if (ged_edit_session_mode(ged_) == GED_EDIT_VIEW) {
      if (!run(std::string(member ? "oed " : "sed ") + quote(path))) return;
    } else {
      if (!ged_edit_session_matches(ged_, path.c_str(), member ? GED_EDIT_INSTANCE : GED_EDIT_PRIMITIVE)) {
        fail("Accept or reject the active edit before transforming another selection");
        return;
      }
    }
    const char* operation = scaling ? "sca " : widget == rotate_ ? "rot " : "tra -m ";
    run(std::string(operation) + value);
  }

  void ged(struct ged* context) { ged_ = context; }

  Fl_Ged_Transform_Editor::CommandHandler command_handler_;
  Fl_Ged_Transform_Editor::ChangedHandler changed_handler_;
  Fl_Ged_Transform_Editor::ErrorHandler error_handler_;
  std::string last_error_;

private:
  Fl_Input* input(const char* label, int row)
  {
    return new Fl_Input(kMargin + kLabelWidth,
                        kMargin + row * (kRowHeight + kGap),
                        kWindowWidth - 2 * kMargin - kLabelWidth,
                        kRowHeight, label);
  }

  bool run(const std::string& command)
  {
    if (!command_handler_) {
      fail("No GED command handler is configured");
      return false;
    }
    if (command_handler_(command) != BRLCAD_OK) {
      fail("The transform failed; see the terminal for details");
      return false;
    }
    last_error_.clear();
    if (changed_handler_) {
      changed_handler_();
    }
    return true;
  }

  void fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
  }

  static void button_callback(Fl_Widget* widget, void* data)
  {
    static_cast<Impl*>(data)->button(widget);
  }

  Fl_Ged_Transform_Editor& owner_;
  struct ged* ged_ = nullptr;
  Fl_Input* path_ = nullptr;
  Fl_Choice* target_ = nullptr;
  Fl_Input* translation_ = nullptr;
  Fl_Input* rotation_ = nullptr;
  Fl_Input* scale_ = nullptr;
  Fl_Button* translate_ = nullptr;
  Fl_Button* rotate_ = nullptr;
  Fl_Button* apply_scale_ = nullptr;
  Fl_Button* reset_ = nullptr;
  Fl_Button* accept_ = nullptr;
  Fl_Button* reject_ = nullptr;
};

Fl_Ged_Transform_Editor::Fl_Ged_Transform_Editor()
  : Fl_Double_Window(kWindowWidth, kWindowHeight, "Object and Matrix Transform"),
    impl_(std::make_unique<Impl>(*this))
{
  set_non_modal();
}
Fl_Ged_Transform_Editor::~Fl_Ged_Transform_Editor() = default;
void Fl_Ged_Transform_Editor::ged(struct ged* context) { impl_->ged(context); }
void Fl_Ged_Transform_Editor::selection(std::string path)
{
  impl_->selection(std::move(path));
}
void Fl_Ged_Transform_Editor::command_handler(CommandHandler handler)
{
  impl_->command_handler_ = std::move(handler);
}
void Fl_Ged_Transform_Editor::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler_ = std::move(handler);
}
void Fl_Ged_Transform_Editor::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}
const std::string& Fl_Ged_Transform_Editor::last_error() const noexcept
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
