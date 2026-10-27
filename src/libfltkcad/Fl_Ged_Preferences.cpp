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

#include "fltkcad/Fl_Ged_Preferences.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Input.H>
#include <FL/fl_ask.H>

#include <array>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace {

constexpr int kWindowWidth = 430;
constexpr int kWindowHeight = 326;
constexpr int kMargin = 12;
constexpr int kRowHeight = 26;
constexpr int kRowGap = 8;
constexpr int kLabelWidth = 170;
constexpr int kButtonWidth = 110;
constexpr int kButtonGap = 8;
constexpr int kInputLimit = 64;
constexpr std::size_t kMaximumHistory = 100000;
constexpr int kMaximumScrollback = 1000000;

bool parse_rgb(const char* text, Fl_Cad_View::Color& color)
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  std::array<int, 3> values{};
  std::string trailing;
  if (!(input >> values[0] >> values[1] >> values[2]) || input >> trailing) {
    return false;
  }
  for (int value : values) {
    if (value < 0 || value > std::numeric_limits<unsigned char>::max()) {
      return false;
    }
  }
  color = {static_cast<unsigned char>(values[0]),
           static_cast<unsigned char>(values[1]),
           static_cast<unsigned char>(values[2])};
  return true;
}

std::string format_rgb(const Fl_Cad_View::Color& color)
{
  return std::to_string(color.red) + " " + std::to_string(color.green) +
         " " + std::to_string(color.blue);
}

template <typename Value>
bool parse_positive(const char* text, Value maximum, Value& value)
{
  if (!text) {
    return false;
  }
  unsigned long long parsed = 0;
  std::istringstream input(text);
  std::string trailing;
  if (!(input >> parsed) || input >> trailing || parsed == 0 ||
      parsed > static_cast<unsigned long long>(maximum)) {
    return false;
  }
  value = static_cast<Value>(parsed);
  return true;
}

} // namespace

class Fl_Ged_Preferences::Impl
{
public:
  explicit Impl(Fl_Ged_Preferences& owner) : owner_(owner)
  {
    owner_.size(kWindowWidth, kWindowHeight);
    owner_.begin();
    background_first_ = make_input("Background top");
    background_second_ = make_input("Background bottom");
    terminal_foreground_ = make_input("Terminal text");
    terminal_background_ = make_input("Terminal background");
    faceplate_font_size_ = make_input("Status font (auto or size)");
    command_history_ = make_input("Command history");
    scrollback_lines_ = make_input("Terminal scrollback");
    wrap_history_ = new Fl_Check_Button(0, 0, 1, kRowHeight, "Wrap historical output");
    cancel_ = new Fl_Button(0, 0, kButtonWidth, kRowHeight, "Cancel");
    apply_ = new Fl_Button(0, 0, kButtonWidth, kRowHeight, "Apply and Save");
    cancel_->callback(cancel_callback, this);
    apply_->callback(apply_callback, this);
    owner_.end();
    owner_.set_modal();
    layout();
    settings(Settings{});
  }

  void settings(const Settings& value)
  {
    settings_ = value;
    background_first_->value(format_rgb(value.background_first).c_str());
    background_second_->value(format_rgb(value.background_second).c_str());
    terminal_foreground_->value(format_rgb(value.terminal_foreground).c_str());
    terminal_background_->value(format_rgb(value.terminal_background).c_str());
    faceplate_font_size_->value(value.faceplate_font_size == 0 ? "auto" :
      std::to_string(value.faceplate_font_size).c_str());
    command_history_->value(std::to_string(value.command_history).c_str());
    scrollback_lines_->value(std::to_string(value.scrollback_lines).c_str());
    wrap_history_->value(value.wrap_history);
    last_error_.clear();
  }

  bool apply()
  {
    Settings replacement;
    if (!parse_rgb(background_first_->value(), replacement.background_first) ||
        !parse_rgb(background_second_->value(),
                   replacement.background_second)) {
      return fail("Background colors must each contain three values from 0 to 255");
    }
    if (!parse_rgb(terminal_foreground_->value(), replacement.terminal_foreground) ||
        !parse_rgb(terminal_background_->value(), replacement.terminal_background)) {
      return fail("Terminal colors must each contain three values from 0 to 255");
    }
    const std::string font = faceplate_font_size_->value();
    if (font != "auto" && font != "0" &&
        !parse_positive(font.c_str(), Fl_Ged_Preferences::maximum_font_size, replacement.faceplate_font_size)) {
      return fail("Status font must be auto or a size between 1 and 256");
    }
    if (!parse_positive(command_history_->value(), kMaximumHistory,
                        replacement.command_history)) {
      return fail("Command history must be between 1 and 100000");
    }
    if (!parse_positive(scrollback_lines_->value(), kMaximumScrollback,
                        replacement.scrollback_lines)) {
      return fail("Terminal scrollback must be between 1 and 1000000 lines");
    }

    replacement.wrap_history = wrap_history_->value() != 0;
    settings_ = replacement;
    last_error_.clear();
    if (changed_handler_) {
      changed_handler_(settings_);
    }
    owner_.hide();
    return true;
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    return false;
  }

  void layout()
  {
    const int input_x = kMargin + kLabelWidth;
    const int input_width = kWindowWidth - input_x - kMargin;
    int row_y = kMargin;
    for (Fl_Input* input : {background_first_, background_second_,
                            terminal_foreground_, terminal_background_,
                            faceplate_font_size_, command_history_, scrollback_lines_}) {
      input->resize(input_x, row_y, input_width, kRowHeight);
      row_y += kRowHeight + kRowGap;
    }
    wrap_history_->resize(kMargin, row_y, kWindowWidth - 2 * kMargin, kRowHeight);
    row_y += kRowHeight + kRowGap;
    const int apply_x = kWindowWidth - kMargin - kButtonWidth;
    apply_->resize(apply_x, row_y, kButtonWidth, kRowHeight);
    cancel_->resize(apply_x - kButtonGap - kButtonWidth, row_y,
                    kButtonWidth, kRowHeight);
  }

  Fl_Input* make_input(const char* label)
  {
    auto* input = new Fl_Input(0, 0, 1, kRowHeight, label);
    input->maximum_size(kInputLimit);
    return input;
  }

  static void apply_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.apply()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void cancel_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->owner_.hide();
  }

  Fl_Ged_Preferences& owner_;
  Fl_Input* background_first_ = nullptr;
  Fl_Input* background_second_ = nullptr;
  Fl_Input* terminal_foreground_ = nullptr;
  Fl_Input* terminal_background_ = nullptr;
  Fl_Input* faceplate_font_size_ = nullptr;
  Fl_Input* command_history_ = nullptr;
  Fl_Input* scrollback_lines_ = nullptr;
  Fl_Check_Button* wrap_history_ = nullptr;
  Fl_Button* cancel_ = nullptr;
  Fl_Button* apply_ = nullptr;
  Settings settings_;
  ChangedHandler changed_handler_;
  std::string last_error_;
};

Fl_Ged_Preferences::Fl_Ged_Preferences(const char* title)
  : Fl_Double_Window(kWindowWidth, kWindowHeight, title),
    impl_(std::make_unique<Impl>(*this))
{
}

Fl_Ged_Preferences::~Fl_Ged_Preferences() = default;

void Fl_Ged_Preferences::settings(const Settings& value)
{
  impl_->settings(value);
}

Fl_Ged_Preferences::Settings Fl_Ged_Preferences::settings() const noexcept
{
  return impl_->settings_;
}

bool Fl_Ged_Preferences::apply() { return impl_->apply(); }

const std::string& Fl_Ged_Preferences::last_error() const noexcept
{
  return impl_->last_error_;
}

void Fl_Ged_Preferences::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler_ = std::move(handler);
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
