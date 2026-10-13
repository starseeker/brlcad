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

#ifndef FLTKCAD_FL_GED_TCL_H
#define FLTKCAD_FL_GED_TCL_H

#include "fltkcad/defines.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

struct ged;

class FLTKCAD_EXPORT Fl_Ged_Tcl
{
public:
  using ChangedHandler = std::function<void()>;
  using OutputHandler = std::function<void(const std::string& text, bool error)>;
  using UiHandler = std::function<int(const std::vector<std::string>& arguments,
                                      std::string& result)>;

  Fl_Ged_Tcl();
  ~Fl_Ged_Tcl();

  Fl_Ged_Tcl(const Fl_Ged_Tcl&) = delete;
  Fl_Ged_Tcl& operator=(const Fl_Ged_Tcl&) = delete;

  bool available() const noexcept;
  // The context is borrowed: unbind before destroying it.  If the application
  // supplies no subprocess handlers, Tcl services the context's process I/O.
  void ged(struct ged* value);
  int evaluate(const std::string& script);
  int source(const std::string& path);
  void interrupt();
  bool command_complete(const std::string& script) const;
  int command_status() const noexcept;
  bool exit_requested() const noexcept;
  int exit_status() const noexcept;
  void output_handler(OutputHandler handler);
  void event_processing(bool enabled);
  bool enable_legacy_tk();
  bool legacy_tk_active() const noexcept;
  const std::string& result() const noexcept;
  const std::string& last_error() const noexcept;
  void changed_handler(ChangedHandler handler);
  void ui_handler(UiHandler handler);

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
