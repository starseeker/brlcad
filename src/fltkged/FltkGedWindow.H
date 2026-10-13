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

#ifndef FLTKGED_FLTK_GED_WINDOW_H
#define FLTKGED_FLTK_GED_WINDOW_H

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Double_Window.H>

#include <memory>
#include <string>

class Fl_Cad_Terminal;
struct ged;

class FltkGedWindow : public Fl_Double_Window
{
public:
  enum class CommandLanguage { Tcl, Ged };
  FltkGedWindow(int width, int height,
                Fl_Cad_View::Backend backend = Fl_Cad_View::Backend::Auto,
                const char* title = "fltkged");
  ~FltkGedWindow() override;

  FltkGedWindow(const FltkGedWindow&) = delete;
  FltkGedWindow& operator=(const FltkGedWindow&) = delete;

  bool open_database(const std::string& path);
  int execute(const std::string& command);
  int execute_prompt(const std::string& command);
  void command_language(CommandLanguage language);
  CommandLanguage command_language() const noexcept;

  struct ged* gedp() noexcept;
  Fl_Cad_View& view() noexcept;
  Fl_Cad_Terminal& terminal() noexcept;
  const std::string& last_error() const noexcept;

  using Fl_Double_Window::show;
  void show() override;
  void hide() override;
  void combined(bool value);
  bool combined() const noexcept;
  Fl_Window& graphics_window() noexcept;

  void resize(int x, int y, int width, int height) override;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

#endif

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
