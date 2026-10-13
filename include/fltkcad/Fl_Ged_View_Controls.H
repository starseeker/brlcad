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

#ifndef FLTKCAD_FL_GED_VIEW_CONTROLS_H
#define FLTKCAD_FL_GED_VIEW_CONTROLS_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <functional>
#include <memory>
#include <string>

class Fl_Cad_View;

class FLTKCAD_EXPORT Fl_Ged_View_Controls : public Fl_Group
{
public:
  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_View_Controls(int x, int y, int width, int height,
                       const char* label = nullptr);
  ~Fl_Ged_View_Controls() override;

  Fl_Ged_View_Controls(const Fl_Ged_View_Controls&) = delete;
  Fl_Ged_View_Controls& operator=(const Fl_Ged_View_Controls&) = delete;

  void view(Fl_Cad_View* view);
  Fl_Cad_View* view() const noexcept;
  bool apply();
  void refresh();

  const std::string& last_error() const noexcept;
  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);
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
