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

#ifndef FLTKCAD_FL_GED_SEARCH_H
#define FLTKCAD_FL_GED_SEARCH_H

#include "fltkcad/defines.h"

#include <FL/Fl_Double_Window.H>

#include <functional>
#include <memory>
#include <string>

struct ged;

class FLTKCAD_EXPORT Fl_Ged_Search : public Fl_Double_Window
{
public:
  using SelectionHandler = std::function<void(const std::string& object_name)>;
  using CommandHandler = std::function<void(const std::string& command)>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_Search();
  ~Fl_Ged_Search() override;

  Fl_Ged_Search(const Fl_Ged_Search&) = delete;
  Fl_Ged_Search& operator=(const Fl_Ged_Search&) = delete;

  void ged(struct ged* value);
  struct ged* ged() const noexcept;
  bool search();
  const std::string& selected_object() const noexcept;
  const std::string& last_error() const noexcept;

  void selection_handler(SelectionHandler handler);
  void command_handler(CommandHandler handler);
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
