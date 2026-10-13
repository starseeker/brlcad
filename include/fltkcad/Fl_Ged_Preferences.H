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

#ifndef FLTKCAD_FL_GED_PREFERENCES_H
#define FLTKCAD_FL_GED_PREFERENCES_H

#include "fltkcad/defines.h"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Double_Window.H>

#include <cstddef>
#include <functional>
#include <memory>
#include <string>

class FLTKCAD_EXPORT Fl_Ged_Preferences : public Fl_Double_Window
{
public:
  static constexpr int maximum_font_size = 256;

  struct Settings
  {
    Fl_Cad_View::Color background_first{};
    Fl_Cad_View::Color background_second{};
    Fl_Cad_View::Color terminal_foreground{};
    Fl_Cad_View::Color terminal_background{255, 255, 255};
    int faceplate_font_size = 0; ///< Zero selects automatic sizing and wrapping.
    std::size_t command_history = 1000;
    int scrollback_lines = 10000;
  };

  using ChangedHandler = std::function<void(const Settings& settings)>;

  explicit Fl_Ged_Preferences(const char* title = "fltkged Preferences");
  ~Fl_Ged_Preferences() override;

  Fl_Ged_Preferences(const Fl_Ged_Preferences&) = delete;
  Fl_Ged_Preferences& operator=(const Fl_Ged_Preferences&) = delete;

  void settings(const Settings& value);
  Settings settings() const noexcept;
  bool apply();

  const std::string& last_error() const noexcept;
  void changed_handler(ChangedHandler handler);

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
