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

#include "common.h"

#include "fltkcad/Fl_Ged_Hierarchy.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Tree.H>

#include <algorithm>
#include <cstring>
#include <utility>

extern "C" {
#include "bu/malloc.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/comb.h"
#include "rt/search.h"
}

namespace {

constexpr int kControlsHeight = 28;
constexpr int kControlGap = 4;
constexpr std::size_t kMaxRootItems = 250000;
constexpr std::size_t kMaxChildrenPerCombination = 250000;
constexpr const char* kLoadingLabel = "Loading...";

std::string boolean_prefix(int operation)
{
  switch (operation) {
    case OP_SUBTRACT: return "- ";
    case OP_INTERSECT: return "+ ";
    case OP_UNION: return "u ";
    default: return {};
  }
}

std::string quoted_argument(const std::string& value)
{
  std::string result;
  result.reserve(value.size() + 2);
  result.push_back('"');
  for (const char character : value) {
    if (character == '\\' || character == '"') {
      result.push_back('\\');
    }
    result.push_back(character);
  }
  result.push_back('"');
  return result;
}

} // namespace

struct Fl_Ged_Hierarchy::Node
{
  std::string object;
  std::vector<std::string> path;
  bool children_loaded = false;
  bool expandable = false;
};

Fl_Ged_Hierarchy::Fl_Ged_Hierarchy(int x, int y, int width, int height,
                                   const char* label)
  : Fl_Group(x, y, width, height, label)
{
  begin();
  tree_ = new Fl_Tree(x, y, width, std::max(1, height - kControlsHeight));
  tree_->showroot(0);
  tree_->selectmode(FL_TREE_SELECT_SINGLE);
  tree_->callback(tree_callback, this);

  const int controls_y = y + std::max(1, height - kControlsHeight) + kControlGap;
  tops_only_ = new Fl_Check_Button(x, controls_y, 92,
                                   kControlsHeight - kControlGap, "Tops only");
  tops_only_->value(1);
  tops_only_->callback(control_callback, this);
  draw_ = new Fl_Button(x + 96, controls_y, 54,
                        kControlsHeight - kControlGap, "Draw");
  draw_->callback(control_callback, this);
  erase_ = new Fl_Button(x + 154, controls_y, 54,
                         kControlsHeight - kControlGap, "Erase");
  erase_->callback(control_callback, this);
  end();
  resizable(tree_);
}

Fl_Ged_Hierarchy::~Fl_Ged_Hierarchy() = default;

void Fl_Ged_Hierarchy::ged(struct ged* value)
{
  ged_ = value;
  refresh();
}

struct ged* Fl_Ged_Hierarchy::ged() const noexcept { return ged_; }

void Fl_Ged_Hierarchy::refresh()
{
  selection_ = {};
  tree_->clear();
  nodes_.clear();
  if (!ged_ || !ged_->dbip) {
    tree_->redraw();
    return;
  }
  populate_roots();
  tree_->redraw();
}

const Fl_Ged_Hierarchy::Selection& Fl_Ged_Hierarchy::selection() const noexcept
{
  return selection_;
}

void Fl_Ged_Hierarchy::selection_handler(SelectionHandler handler)
{
  selection_handler_ = std::move(handler);
}

void Fl_Ged_Hierarchy::command_handler(CommandHandler handler)
{
  command_handler_ = std::move(handler);
}

void Fl_Ged_Hierarchy::error_handler(ErrorHandler handler)
{
  error_handler_ = std::move(handler);
}

void Fl_Ged_Hierarchy::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  const int tree_height = std::max(1, height - kControlsHeight);
  tree_->resize(x, y, width, tree_height);
  const int controls_y = y + tree_height + kControlGap;
  tops_only_->resize(x, controls_y, 92, kControlsHeight - kControlGap);
  draw_->resize(x + 96, controls_y, 54, kControlsHeight - kControlGap);
  erase_->resize(x + 154, controls_y, 54, kControlsHeight - kControlGap);
}

void Fl_Ged_Hierarchy::populate_roots()
{
  db_update_nref(ged_->dbip);
  std::vector<std::string> names;
  names.reserve(std::min<std::size_t>(kMaxRootItems, 1024));
  bool exceeded = false;
  struct directory* entry = nullptr;
  FOR_ALL_DIRECTORY_START(entry, ged_->dbip)
    if (!entry->d_namep || entry->d_addr == RT_DIR_PHONY_ADDR ||
        (entry->d_flags & RT_DIR_HIDDEN) ||
        (tops_only_->value() && entry->d_nref != 0)) {
      continue;
    }
    if (names.size() == kMaxRootItems) {
      exceeded = true;
      continue;
    }
    names.emplace_back(entry->d_namep);
  FOR_ALL_DIRECTORY_END;
  if (exceeded) {
    report_error("Database hierarchy exceeds the configured root item limit");
    return;
  }
  std::sort(names.begin(), names.end());

  Fl_Tree_Item* root = tree_->root();
  if (!root) {
    return;
  }
  for (const auto& name : names) {
    add_node(*root, name, {name}, OP_UNION);
  }
}

void Fl_Ged_Hierarchy::populate_children(Fl_Tree_Item& item, Node& node)
{
  if (node.children_loaded || !node.expandable || !ged_ || !ged_->dbip) {
    return;
  }
  node.children_loaded = true;
  item.clear_children();

  struct directory* parent =
    db_lookup(ged_->dbip, node.object.c_str(), LOOKUP_QUIET);
  if (!parent || !(parent->d_flags & RT_DIR_COMB)) {
    return;
  }

  struct rt_db_internal internal;
  RT_DB_INTERNAL_INIT(&internal);
  if (rt_db_get_internal(&internal, parent, ged_->dbip, nullptr) < 0) {
    report_error("Unable to read combination " + node.object);
    return;
  }
  auto* combination = static_cast<struct rt_comb_internal*>(internal.idb_ptr);
  struct directory** children = nullptr;
  int* operations = nullptr;
  const int count =
    db_comb_children(ged_->dbip, combination, &children, &operations, nullptr);
  rt_db_free_internal(&internal);
  if (count < 0 || static_cast<std::size_t>(count) >
                     kMaxChildrenPerCombination) {
    if (children) {
      bu_free(children, "combination children");
    }
    if (operations) {
      bu_free(operations, "combination operations");
    }
    report_error("Combination child count exceeds the configured limit");
    return;
  }

  for (int index = 0; index < count; ++index) {
    if (!children[index] || !children[index]->d_namep) {
      continue;
    }
    const std::string child(children[index]->d_namep);
    auto path = node.path;
    path.push_back(child);
    const int operation = operations ? operations[index] : OP_UNION;
    add_node(item, child, std::move(path), operation);
    if (is_cycle(node.path, child)) {
      auto* added = static_cast<Node*>(item.child(item.children() - 1)->user_data());
      if (added) {
        added->expandable = false;
      }
      item.child(item.children() - 1)->clear_children();
    }
  }
  if (children) {
    bu_free(children, "combination children");
  }
  if (operations) {
    bu_free(operations, "combination operations");
  }
}

void Fl_Ged_Hierarchy::add_node(Fl_Tree_Item& parent,
                                const std::string& object,
                                std::vector<std::string> path,
                                int boolean_operation)
{
  auto node = std::make_unique<Node>();
  node->object = object;
  node->path = std::move(path);
  if (ged_ && ged_->dbip) {
    struct directory* entry =
      db_lookup(ged_->dbip, object.c_str(), LOOKUP_QUIET);
    node->expandable = entry && (entry->d_flags & RT_DIR_COMB);
  }

  const std::string text = boolean_prefix(boolean_operation) + object;
  Fl_Tree_Item* item = parent.add(tree_->prefs(), text.c_str());
  if (!item) {
    report_error("Unable to add an item to the database hierarchy");
    return;
  }
  item->user_data(node.get());
  if (node->expandable) {
    item->add(tree_->prefs(), kLoadingLabel);
    item->close();
  }
  nodes_.push_back(std::move(node));
}

void Fl_Ged_Hierarchy::update_selection(Fl_Tree_Item* item)
{
  auto* node = item ? static_cast<Node*>(item->user_data()) : nullptr;
  if (!node) {
    return;
  }
  selection_.object = node->object;
  selection_.path = node->path;
  if (selection_handler_) {
    selection_handler_(selection_);
  }
}

void Fl_Ged_Hierarchy::report_error(std::string message)
{
  if (error_handler_) {
    error_handler_(message);
  }
}

bool Fl_Ged_Hierarchy::is_cycle(const std::vector<std::string>& path,
                                const std::string& child) const
{
  return std::find(path.begin(), path.end(), child) != path.end();
}

void Fl_Ged_Hierarchy::tree_callback(Fl_Widget*, void* data)
{
  auto& self = *static_cast<Fl_Ged_Hierarchy*>(data);
  Fl_Tree_Item* item = self.tree_->callback_item();
  if (!item) {
    return;
  }
  if (self.tree_->callback_reason() == FL_TREE_REASON_OPENED) {
    auto* node = static_cast<Node*>(item->user_data());
    if (node) {
      self.populate_children(*item, *node);
    }
  }
  if (self.tree_->callback_reason() == FL_TREE_REASON_SELECTED) {
    self.update_selection(item);
  }
}

void Fl_Ged_Hierarchy::control_callback(Fl_Widget* widget, void* data)
{
  auto& self = *static_cast<Fl_Ged_Hierarchy*>(data);
  if (widget == self.tops_only_) {
    self.refresh();
    return;
  }
  if (!self.command_handler_ || self.selection_.object.empty()) {
    return;
  }
  if (widget == self.draw_) {
    self.command_handler_("draw " + quoted_argument(self.selection_.object));
  } else if (widget == self.erase_) {
    self.command_handler_("erase " + quoted_argument(self.selection_.object));
  }
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
