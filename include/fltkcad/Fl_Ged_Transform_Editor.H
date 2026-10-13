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

#ifndef FLTKCAD_FL_GED_TRANSFORM_EDITOR_H
#define FLTKCAD_FL_GED_TRANSFORM_EDITOR_H

#include "fltkcad/defines.h"

#include <FL/Fl_Double_Window.H>

#include <functional>
#include <memory>
#include <string>

struct ged;

class FLTKCAD_EXPORT Fl_Ged_Transform_Editor : public Fl_Double_Window
{
public:
  using CommandHandler = std::function<int(const std::string& command)>;
  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_Transform_Editor();
  ~Fl_Ged_Transform_Editor() override;

  Fl_Ged_Transform_Editor(const Fl_Ged_Transform_Editor&) = delete;
  Fl_Ged_Transform_Editor& operator=(const Fl_Ged_Transform_Editor&) = delete;

  void ged(struct ged* context);
  void selection(std::string path);
  void command_handler(CommandHandler handler);
  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);
  const std::string& last_error() const noexcept;

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
