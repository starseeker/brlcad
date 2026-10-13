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

#ifndef FLTKCAD_FL_GED_HIERARCHY_H
#define FLTKCAD_FL_GED_HIERARCHY_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <functional>
#include <memory>
#include <string>
#include <vector>

class Fl_Button;
class Fl_Check_Button;
class Fl_Tree;
class Fl_Tree_Item;
struct ged;

class FLTKCAD_EXPORT Fl_Ged_Hierarchy : public Fl_Group
{
public:
  struct Selection
  {
    std::string object;
    std::vector<std::string> path;
  };

  using SelectionHandler = std::function<void(const Selection& selection)>;
  using CommandHandler = std::function<void(const std::string& command)>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Ged_Hierarchy(int x, int y, int width, int height,
                   const char* label = nullptr);
  ~Fl_Ged_Hierarchy() override;

  Fl_Ged_Hierarchy(const Fl_Ged_Hierarchy&) = delete;
  Fl_Ged_Hierarchy& operator=(const Fl_Ged_Hierarchy&) = delete;

  void ged(struct ged* value);
  struct ged* ged() const noexcept;
  void refresh();

  const Selection& selection() const noexcept;
  void selection_handler(SelectionHandler handler);
  void command_handler(CommandHandler handler);
  void error_handler(ErrorHandler handler);

  void resize(int x, int y, int width, int height) override;

private:
  struct Node;
  void populate_roots();
  void populate_children(Fl_Tree_Item& item, Node& node);
  void add_node(Fl_Tree_Item& parent, const std::string& object,
                std::vector<std::string> path, int boolean_operation);
  void update_selection(Fl_Tree_Item* item);
  void report_error(std::string message);
  bool is_cycle(const std::vector<std::string>& path,
                const std::string& child) const;
  static void tree_callback(Fl_Widget*, void* data);
  static void control_callback(Fl_Widget*, void* data);

  struct ged* ged_ = nullptr;
  Fl_Tree* tree_ = nullptr;
  Fl_Check_Button* tops_only_ = nullptr;
  Fl_Button* draw_ = nullptr;
  Fl_Button* erase_ = nullptr;
  std::vector<std::unique_ptr<Node>> nodes_;
  Selection selection_;
  SelectionHandler selection_handler_;
  CommandHandler command_handler_;
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
