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

#ifndef FLTKCAD_FL_GED_ATTRIBUTES_H
#define FLTKCAD_FL_GED_ATTRIBUTES_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <functional>
#include <string>
#include <utility>
#include <vector>

class Fl_Button;
class Fl_Table_Row;
struct ged;

class FLTKCAD_EXPORT Fl_Ged_Attributes : public Fl_Group
{
public:
  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_Attributes(int x, int y, int width, int height,
                    const char* label = nullptr);
  ~Fl_Ged_Attributes() override;

  Fl_Ged_Attributes(const Fl_Ged_Attributes&) = delete;
  Fl_Ged_Attributes& operator=(const Fl_Ged_Attributes&) = delete;

  void ged(struct ged* value);
  void object(std::string name);
  const std::string& object() const noexcept;
  void refresh();

  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);

  void resize(int x, int y, int width, int height) override;

private:
  class Table;
  bool write_attributes();
  void add_attribute();
  void edit_attribute();
  void remove_attribute();
  void update_controls();
  void report_error(std::string message);
  static void button_callback(Fl_Widget* widget, void* data);

  struct ged* ged_ = nullptr;
  std::string object_;
  std::vector<std::pair<std::string, std::string>> attributes_;
  bool attributes_complete_ = true;
  Table* table_ = nullptr;
  Fl_Button* add_ = nullptr;
  Fl_Button* edit_ = nullptr;
  Fl_Button* remove_ = nullptr;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
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
