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

#ifndef FLTKCAD_FL_GED_INTERACTION_H
#define FLTKCAD_FL_GED_INTERACTION_H

#include "fltkcad/defines.h"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Group.H>

#include <functional>
#include <memory>
#include <string>

struct ged;

class FLTKCAD_EXPORT Fl_Ged_Interaction : public Fl_Group
{
public:
  enum class Mode
  {
    Navigate,
    Select,
    QueryRay,
    RectangleSelect,
    RectangleZoom,
    TranslateObject,
    RotateObject,
    ScaleObject
  };

  struct Selection
  {
    std::string path;
    std::string object;
  };

  using CommandHandler = std::function<int(const std::string& command)>;
  using SelectionHandler = std::function<void(const Selection& selection)>;
  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_Interaction(int x, int y, int width, int height,
                     const char* label = nullptr);
  ~Fl_Ged_Interaction() override;

  Fl_Ged_Interaction(const Fl_Ged_Interaction&) = delete;
  Fl_Ged_Interaction& operator=(const Fl_Ged_Interaction&) = delete;

  void ged(struct ged* value);
  void selection(std::string path, std::string object);
  Mode mode() const noexcept;
  void mode(Mode value);

  bool process_input(Fl_Cad_View& view,
                     const Fl_Cad_View::InputEvent& event);
  const std::string& last_error() const noexcept;

  void command_handler(CommandHandler handler);
  void selection_handler(SelectionHandler handler);
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
