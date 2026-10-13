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

#include "fltkcad/Fl_Ged_Render_Tools.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Native_File_Chooser.H>

#include <algorithm>
#include <string>
#include <utility>

extern "C" {
#include "bu/defines.h"
}

namespace {

constexpr int kWindowWidth = 520;
constexpr int kWindowHeight = 246;
constexpr int kMargin = 12;
constexpr int kRowHeight = 28;
constexpr int kGap = 8;
constexpr int kLabelWidth = 108;
constexpr int kOptionsLimit = 4096;
constexpr int kColorLimit = 64;

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

} // namespace

class Fl_Ged_Render_Tools::Impl
{
public:
  explicit Impl(Fl_Ged_Render_Tools& owner) : owner_(owner)
  {
    owner_.begin();
    options_ = new Fl_Input(kMargin + kLabelWidth, kMargin,
                            kWindowWidth - 2 * kMargin - kLabelWidth,
                            kRowHeight, "RT options");
    options_->value("-P 1");
    options_->maximum_size(kOptionsLimit);

    const int render_y = kMargin + kRowHeight + kGap;
    embedded_ = new Fl_Button(kMargin, render_y, 155, kRowHeight,
                              "Render in View");
    embedded_->callback(button_callback, this);
    file_ = new Fl_Button(kMargin + 163, render_y, 155, kRowHeight,
                          "Render to PIX...");
    file_->callback(button_callback, this);
    abort_ = new Fl_Button(kMargin + 326, render_y, 155, kRowHeight,
                           "Abort Render");
    abort_->callback(button_callback, this);

    const int query_y = render_y + kRowHeight + kGap;
    center_ray_ = new Fl_Button(kMargin, query_y, 236, kRowHeight,
                                "Query Center Ray");
    center_ray_->callback(button_callback, this);
    clear_ = new Fl_Button(kMargin + 244, query_y, 237, kRowHeight,
                           "Clear Framebuffer");
    clear_->callback(button_callback, this);

    const int framebuffer_y = query_y + kRowHeight + kGap;
    load_ = new Fl_Button(kMargin, framebuffer_y, 236, kRowHeight,
                          "Load PNG into Framebuffer...");
    load_->callback(button_callback, this);
    save_ = new Fl_Button(kMargin + 244, framebuffer_y, 237, kRowHeight,
                          "Save Framebuffer as PIX...");
    save_->callback(button_callback, this);

    const int color_y = framebuffer_y + kRowHeight + kGap;
    clear_color_ = new Fl_Input(kMargin + kLabelWidth, color_y,
                                kWindowWidth - 2 * kMargin - kLabelWidth,
                                kRowHeight, "Clear RGB");
    clear_color_->value("0 0 0");
    clear_color_->maximum_size(kColorLimit);
    owner_.end();
    owner_.resizable(options_);
  }

  void run(const std::string& command)
  {
    if (!command_handler_) {
      fail("No GED command handler is configured");
      return;
    }
    if (command_handler_(command) != BRLCAD_OK) {
      fail("The requested render operation failed; see the terminal for details");
      return;
    }
    last_error_.clear();
  }

  void button(Fl_Widget* widget)
  {
    const std::string options = options_->value() ? options_->value() : "";
    if (widget == embedded_) {
      run(options.empty() ? "ert" : "ert " + options);
      return;
    }
    if (widget == abort_) {
      run("rtabort");
      return;
    }
    if (widget == center_ray_) {
      run("vquery_ray 0 0");
      return;
    }
    if (widget == clear_) {
      const std::string color =
        clear_color_->value() ? clear_color_->value() : "";
      run(color.empty() ? "fbclear" : "fbclear " + quote(color));
      return;
    }
    if (widget == file_) {
      choose_and_run("Raytrace to a PIX image", "PIX image\t*.pix", "rt",
                     options, " -o ", Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
      return;
    }
    if (widget == load_) {
      choose_and_run("Load a PNG image", "PNG image\t*.png", "png2fb", {},
                     " ", Fl_Native_File_Chooser::BROWSE_FILE);
      return;
    }
    if (widget == save_) {
      choose_and_run("Save framebuffer as PIX", "PIX image\t*.pix", "fb2pix",
                     {}, " ", Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
    }
  }

  Fl_Ged_Render_Tools::CommandHandler command_handler_;
  Fl_Ged_Render_Tools::ErrorHandler error_handler_;
  std::string last_error_;

private:
  void choose_and_run(const char* title, const char* filter,
                      const std::string& command, const std::string& options,
                      const char* separator,
                      Fl_Native_File_Chooser::Type type)
  {
    Fl_Native_File_Chooser chooser;
    chooser.title(title);
    chooser.type(type);
    chooser.filter(filter);
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }
    std::string invocation = command;
    if (!options.empty()) {
      invocation += " " + options;
    }
    invocation += separator;
    invocation += quote(chooser.filename());
    run(invocation);
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

  Fl_Ged_Render_Tools& owner_;
  Fl_Input* options_ = nullptr;
  Fl_Input* clear_color_ = nullptr;
  Fl_Button* embedded_ = nullptr;
  Fl_Button* file_ = nullptr;
  Fl_Button* abort_ = nullptr;
  Fl_Button* center_ray_ = nullptr;
  Fl_Button* clear_ = nullptr;
  Fl_Button* load_ = nullptr;
  Fl_Button* save_ = nullptr;
};

Fl_Ged_Render_Tools::Fl_Ged_Render_Tools()
  : Fl_Double_Window(kWindowWidth, kWindowHeight, "Raytrace and Framebuffer"),
    impl_(std::make_unique<Impl>(*this))
{
  set_non_modal();
}

Fl_Ged_Render_Tools::~Fl_Ged_Render_Tools() = default;
void Fl_Ged_Render_Tools::command_handler(CommandHandler handler)
{
  impl_->command_handler_ = std::move(handler);
}
void Fl_Ged_Render_Tools::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}
const std::string& Fl_Ged_Render_Tools::last_error() const noexcept
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
