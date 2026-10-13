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

#ifndef FLTKCAD_FL_GED_HELP_H
#define FLTKCAD_FL_GED_HELP_H

#include "fltkcad/defines.h"

#include <FL/Fl_Double_Window.H>

#include <memory>
#include <string>

class Fl_Litehtml_View;

class FLTKCAD_EXPORT Fl_Ged_Help : public Fl_Double_Window
{
public:
  Fl_Ged_Help(int width = 900, int height = 700,
              const char* title = "MGED Documentation");
  ~Fl_Ged_Help() override;

  Fl_Ged_Help(const Fl_Ged_Help&) = delete;
  Fl_Ged_Help& operator=(const Fl_Ged_Help&) = delete;

  bool load_default();
  bool manuals(const std::string& language = "en", char section = '\0');
  std::size_t manual_count() const noexcept;
  bool back();
  bool forward();
  bool load(const std::string& path);
  const std::string& last_error() const noexcept;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

class FLTKCAD_EXPORT Fl_Ged_About : public Fl_Double_Window
{
public:
  explicit Fl_Ged_About(const char* backend_name = nullptr);
  ~Fl_Ged_About() override;

private:
  static void close_callback(Fl_Widget*, void* data);
  static void docs_callback(Fl_Widget*, void* data);

  std::unique_ptr<Fl_Ged_Help> help_;
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
