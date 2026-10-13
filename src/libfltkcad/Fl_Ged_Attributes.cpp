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

#include "fltkcad/Fl_Ged_Attributes.H"

#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Table_Row.H>
#include <FL/fl_ask.H>
#include <FL/fl_draw.H>

#include <algorithm>
#include <cstring>
#include <utility>

extern "C" {
#include "bu/avs.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/db_attr.h"
}

namespace {

constexpr int kControlsHeight = 30;
constexpr int kButtonGap = 4;
constexpr int kButtonWidth = 64;
constexpr int kHeaderHeight = 24;
constexpr int kRowHeight = 22;
constexpr std::size_t kMaxAttributes = 100000;
constexpr std::size_t kMaxAttributeNameBytes = 4096;
constexpr std::size_t kMaxAttributeValueBytes = 64U * 1024U;

} // namespace

class Fl_Ged_Attributes::Table : public Fl_Table_Row
{
public:
  Table(int x, int y, int width, int height,
        const std::vector<std::pair<std::string, std::string>>& values)
    : Fl_Table_Row(x, y, width, height), values_(values)
  {
    cols(2);
    col_header(1);
    col_resize(1);
    col_width(0, std::max(80, width / 3));
    col_width(1, std::max(100, width - col_width(0) - 4));
    row_header(0);
    row_resize(0);
    row_height_all(kRowHeight);
    col_header_height(kHeaderHeight);
    type(SELECT_SINGLE);
    end();
  }

  void sync()
  {
    rows(static_cast<int>(values_.size()));
    redraw();
  }

  int selected_row()
  {
    for (int row = 0; row < rows(); ++row) {
      if (row_selected(row)) {
        return row;
      }
    }
    return -1;
  }

protected:
  void draw_cell(TableContext context, int row, int column, int x, int y,
                 int width, int height) override
  {
    switch (context) {
      case CONTEXT_COL_HEADER:
        fl_push_clip(x, y, width, height);
        fl_draw_box(FL_THIN_UP_BOX, x, y, width, height, FL_BACKGROUND_COLOR);
        fl_color(FL_FOREGROUND_COLOR);
        fl_draw(column == 0 ? "Attribute" : "Value", x + 4, y, width - 8,
                height, FL_ALIGN_LEFT);
        fl_pop_clip();
        break;
      case CONTEXT_CELL: {
        fl_push_clip(x, y, width, height);
        const Fl_Color background = row_selected(row) ? selection_color()
                                                       : FL_BACKGROUND2_COLOR;
        fl_color(background);
        fl_rectf(x, y, width, height);
        fl_color(FL_LIGHT2);
        fl_rect(x, y, width, height);
        fl_color(fl_contrast(FL_FOREGROUND_COLOR, background));
        if (row >= 0 && static_cast<std::size_t>(row) < values_.size()) {
          const auto& value = values_[static_cast<std::size_t>(row)];
          const std::string& text = column == 0 ? value.first : value.second;
          fl_draw(text.c_str(), x + 4, y, width - 8, height, FL_ALIGN_LEFT);
        }
        fl_pop_clip();
        break;
      }
      default: break;
    }
  }

private:
  const std::vector<std::pair<std::string, std::string>>& values_;
};

Fl_Ged_Attributes::Fl_Ged_Attributes(int x, int y, int width, int height,
                                     const char* label)
  : Fl_Group(x, y, width, height, label)
{
  begin();
  table_ = new Table(x, y, width, std::max(1, height - kControlsHeight),
                     attributes_);
  const int controls_y = y + std::max(1, height - kControlsHeight) + kButtonGap;
  add_ = new Fl_Button(x, controls_y, kButtonWidth,
                       kControlsHeight - kButtonGap, "Add");
  edit_ = new Fl_Button(x + kButtonWidth + kButtonGap, controls_y,
                        kButtonWidth, kControlsHeight - kButtonGap, "Edit");
  remove_ = new Fl_Button(x + (kButtonWidth + kButtonGap) * 2, controls_y,
                          kButtonWidth, kControlsHeight - kButtonGap, "Delete");
  add_->callback(button_callback, this);
  edit_->callback(button_callback, this);
  remove_->callback(button_callback, this);
  end();
  resizable(table_);
  update_controls();
}

Fl_Ged_Attributes::~Fl_Ged_Attributes() = default;

void Fl_Ged_Attributes::ged(struct ged* value)
{
  ged_ = value;
  refresh();
}

void Fl_Ged_Attributes::object(std::string name)
{
  if (object_ == name) {
    return;
  }
  object_ = std::move(name);
  refresh();
}

const std::string& Fl_Ged_Attributes::object() const noexcept
{
  return object_;
}

void Fl_Ged_Attributes::refresh()
{
  attributes_.clear();
  attributes_complete_ = true;
  if (!ged_ || !ged_->dbip || object_.empty()) {
    table_->sync();
    update_controls();
    return;
  }

  struct directory* entry =
    db_lookup(ged_->dbip, object_.c_str(), LOOKUP_QUIET);
  if (!entry) {
    object_.clear();
    table_->sync();
    update_controls();
    report_error("The selected database object no longer exists");
    return;
  }

  struct bu_attribute_value_set values = BU_AVS_INIT_ZERO;
  bu_avs_init_empty(&values);
  if (db5_get_attributes(ged_->dbip, &values, entry) < 0) {
    attributes_complete_ = false;
    bu_avs_free(&values);
    table_->sync();
    update_controls();
    report_error("Unable to read attributes for " + object_);
    return;
  }

  struct bu_attribute_value_pair* pair = nullptr;
  for (BU_AVS_FOR(pair, &values)) {
    if (attributes_.size() >= kMaxAttributes) {
      attributes_complete_ = false;
      break;
    }
    const char* name = pair->name ? pair->name : "";
    const char* value = pair->value ? pair->value : "";
    const std::size_t name_length = std::strlen(name);
    const std::size_t value_length = std::strlen(value);
    if (name_length > kMaxAttributeNameBytes ||
        value_length > kMaxAttributeValueBytes) {
      attributes_complete_ = false;
      continue;
    }
    attributes_.emplace_back(name, value);
  }
  bu_avs_free(&values);
  std::sort(attributes_.begin(), attributes_.end(),
            [](const auto& left, const auto& right) {
              return left.first < right.first;
            });
  table_->sync();
  update_controls();
  if (!attributes_complete_) {
    report_error(
      "Some attributes exceed the configured display limits; editing is disabled");
  }
}

void Fl_Ged_Attributes::changed_handler(ChangedHandler handler)
{
  changed_handler_ = std::move(handler);
}

void Fl_Ged_Attributes::error_handler(ErrorHandler handler)
{
  error_handler_ = std::move(handler);
}

void Fl_Ged_Attributes::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  const int table_height = std::max(1, height - kControlsHeight);
  table_->resize(x, y, width, table_height);
  table_->col_width(0, std::max(80, width / 3));
  table_->col_width(1, std::max(100, width - table_->col_width(0) - 4));
  const int controls_y = y + table_height + kButtonGap;
  add_->resize(x, controls_y, kButtonWidth, kControlsHeight - kButtonGap);
  edit_->resize(x + kButtonWidth + kButtonGap, controls_y, kButtonWidth,
                kControlsHeight - kButtonGap);
  remove_->resize(x + (kButtonWidth + kButtonGap) * 2, controls_y,
                  kButtonWidth, kControlsHeight - kButtonGap);
}

bool Fl_Ged_Attributes::write_attributes()
{
  if (!ged_ || !ged_->dbip || object_.empty()) {
    report_error("No database object is selected");
    return false;
  }
  if (!attributes_complete_) {
    report_error(
      "Attributes cannot be written while the displayed set is incomplete");
    return false;
  }
  if (ged_->dbip->dbi_read_only) {
    report_error("The database is read-only");
    return false;
  }
  if (db_version(ged_->dbip) < 5) {
    report_error("Editing attributes requires a version 5 database");
    return false;
  }
  struct directory* entry =
    db_lookup(ged_->dbip, object_.c_str(), LOOKUP_QUIET);
  if (!entry) {
    report_error("The selected database object no longer exists");
    return false;
  }

  struct bu_attribute_value_set values = BU_AVS_INIT_ZERO;
  bu_avs_init(&values, attributes_.size(), "FLTK attribute editor");
  for (const auto& pair : attributes_) {
    bu_avs_add(&values, pair.first.c_str(), pair.second.c_str());
  }
  if (db5_replace_attributes(entry, &values, ged_->dbip) < 0) {
    report_error("Unable to write attributes for " + object_);
    return false;
  }
  if (changed_handler_) {
    changed_handler_();
  }
  refresh();
  return true;
}

void Fl_Ged_Attributes::add_attribute()
{
  const char* name = fl_input("Attribute name:");
  if (!name) {
    return;
  }
  const std::string attribute_name(name);
  if (attribute_name.empty() || attribute_name.size() > kMaxAttributeNameBytes ||
      attribute_name.find('\0') != std::string::npos) {
    report_error("Attribute name is empty or exceeds the configured limit");
    return;
  }
  const auto duplicate = std::find_if(
    attributes_.begin(), attributes_.end(), [&](const auto& pair) {
      return pair.first == attribute_name;
    });
  if (duplicate != attributes_.end()) {
    report_error("An attribute with that name already exists");
    return;
  }
  const char* value = fl_input("Attribute value:");
  if (!value) {
    return;
  }
  if (std::strlen(value) > kMaxAttributeValueBytes) {
    report_error("Attribute value exceeds the configured limit");
    return;
  }
  attributes_.emplace_back(attribute_name, value);
  if (!write_attributes()) {
    attributes_.pop_back();
  }
}

void Fl_Ged_Attributes::edit_attribute()
{
  const int row = table_->selected_row();
  if (row < 0 || static_cast<std::size_t>(row) >= attributes_.size()) {
    report_error("Select an attribute to edit");
    return;
  }
  auto& attribute = attributes_[static_cast<std::size_t>(row)];
  const std::string prompt = "Value for " + attribute.first + ':';
  const char* value = fl_input("%s", attribute.second.c_str(), prompt.c_str());
  if (!value) {
    return;
  }
  if (std::strlen(value) > kMaxAttributeValueBytes) {
    report_error("Attribute value exceeds the configured limit");
    return;
  }
  const std::string old_value = attribute.second;
  attribute.second = value;
  if (!write_attributes()) {
    attribute.second = old_value;
  }
}

void Fl_Ged_Attributes::remove_attribute()
{
  const int row = table_->selected_row();
  if (row < 0 || static_cast<std::size_t>(row) >= attributes_.size()) {
    report_error("Select an attribute to delete");
    return;
  }
  const auto index = static_cast<std::size_t>(row);
  const std::string prompt =
    "Delete attribute " + attributes_[index].first + '?';
  if (fl_choice("%s", "Cancel", "Delete", nullptr, prompt.c_str()) != 1) {
    return;
  }
  const auto removed = attributes_[index];
  attributes_.erase(attributes_.begin() + row);
  if (!write_attributes()) {
    attributes_.insert(attributes_.begin() + row, removed);
  }
}

void Fl_Ged_Attributes::update_controls()
{
  const bool writable = ged_ && ged_->dbip && !object_.empty() &&
                        attributes_complete_ && !ged_->dbip->dbi_read_only &&
                        db_version(ged_->dbip) >= 5;
  if (writable) {
    add_->activate();
    edit_->activate();
    remove_->activate();
  } else {
    add_->deactivate();
    edit_->deactivate();
    remove_->deactivate();
  }
}

void Fl_Ged_Attributes::report_error(std::string message)
{
  if (error_handler_) {
    error_handler_(message);
  }
}

void Fl_Ged_Attributes::button_callback(Fl_Widget* widget, void* data)
{
  auto& self = *static_cast<Fl_Ged_Attributes*>(data);
  if (widget == self.add_) {
    self.add_attribute();
  } else if (widget == self.edit_) {
    self.edit_attribute();
  } else if (widget == self.remove_) {
    self.remove_attribute();
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
