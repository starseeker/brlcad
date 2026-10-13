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

#include "fltkcad/Fl_Ged_Combination_Editor.H"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Multiline_Input.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bu/vls.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/geom.h"
}

namespace {

constexpr int kInitialWidth = 780;
constexpr int kInitialHeight = 600;
constexpr int kMargin = 10;
constexpr int kRowHeight = 28;
constexpr int kLabelWidth = 105;
constexpr int kColumnGap = 12;
constexpr int kButtonWidth = 90;
constexpr int kButtonHeight = 30;
constexpr std::size_t kMaxTreeBytes = 1024U * 1024U;
constexpr std::size_t kMaxTextBytes = 64U * 1024U;
constexpr long kMaxColorChannel = 255;

bool parse_integer(const char* text, long& value)
{
  if (!text || !*text) {
    return false;
  }
  errno = 0;
  char* end = nullptr;
  const long parsed = std::strtol(text, &end, 10);
  if (errno == ERANGE || text == end || *end != '\0') {
    return false;
  }
  value = parsed;
  return true;
}

std::string integer_text(long value)
{
  return std::to_string(value);
}

bool valid_color(const char* text)
{
  if (!text || !*text) {
    return true;
  }
  std::istringstream input(text);
  long red = 0;
  long green = 0;
  long blue = 0;
  std::string trailing;
  if (!(input >> red >> green >> blue) || (input >> trailing)) {
    return false;
  }
  return red >= 0 && red <= kMaxColorChannel && green >= 0 &&
         green <= kMaxColorChannel && blue >= 0 && blue <= kMaxColorChannel;
}

} // namespace

class Fl_Ged_Combination_Editor::Impl
{
public:
  explicit Impl(Fl_Ged_Combination_Editor& owner) : owner_(owner)
  {
    owner_.begin();
    heading_ = new Fl_Box(kMargin, kMargin, owner_.w() - 2 * kMargin,
                          kRowHeight, "No combination loaded");
    heading_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);

    region_ = new Fl_Check_Button(kMargin + kLabelWidth,
                                  kMargin + kRowHeight, 110, kRowHeight,
                                  "Region");
    inherit_ = new Fl_Check_Button(kMargin + kLabelWidth + 125,
                                   kMargin + kRowHeight, 110, kRowHeight,
                                   "Inherit");

    region_id_ = make_input("Region ID", 2);
    air_ = make_input("Air code", 3);
    los_ = make_input("LOS", 4);
    material_id_ = make_input("Material ID", 5);
    color_ = make_input("Color RGB", 6);
    shader_ = make_input("Shader", 7);
    material_ = make_input("Material", 8);

    tree_ = new Fl_Multiline_Input(kMargin + kLabelWidth,
                                   kMargin + 9 * kRowHeight,
                                   owner_.w() - kLabelWidth - 2 * kMargin,
                                   owner_.h() - 10 * kRowHeight -
                                     kButtonHeight - 3 * kMargin,
                                   "Boolean tree");
    tree_->align(FL_ALIGN_LEFT | FL_ALIGN_TOP);
    add_union_ = new Fl_Button(kMargin,
                               owner_.h() - kButtonHeight - kMargin,
                               kButtonWidth, kButtonHeight, "Add Union");
    add_subtract_ = new Fl_Button(kMargin + kButtonWidth + kColumnGap,
                                  owner_.h() - kButtonHeight - kMargin,
                                  kButtonWidth, kButtonHeight, "Add Subtract");
    add_intersect_ = new Fl_Button(kMargin + 2 * (kButtonWidth + kColumnGap),
                                   owner_.h() - kButtonHeight - kMargin,
                                   kButtonWidth, kButtonHeight, "Add Intersect");
    save_ = new Fl_Button(owner_.w() - 3 * (kButtonWidth + kColumnGap),
                          owner_.h() - kButtonHeight - kMargin, kButtonWidth,
                          kButtonHeight, "Save");
    reload_ = new Fl_Button(owner_.w() - 2 * (kButtonWidth + kColumnGap),
                            owner_.h() - kButtonHeight - kMargin, kButtonWidth,
                            kButtonHeight, "Reload");
    close_ = new Fl_Button(owner_.w() - kButtonWidth - kMargin,
                           owner_.h() - kButtonHeight - kMargin, kButtonWidth,
                           kButtonHeight, "Close");
    save_->callback(save_callback, this);
    reload_->callback(reload_callback, this);
    close_->callback(close_callback, this);
    add_union_->callback(add_member_callback, this);
    add_subtract_->callback(add_member_callback, this);
    add_intersect_->callback(add_member_callback, this);
    region_->callback(modified_callback, this);
    inherit_->callback(modified_callback, this);
    Fl_Input* inputs[] = {region_id_, air_, los_, material_id_, color_,
                          shader_, material_, tree_};
    for (auto* input : inputs) {
      input->when(FL_WHEN_CHANGED);
      input->callback(modified_callback, this);
    }
    owner_.callback(close_callback, this);
    owner_.end();
    owner_.resizable(tree_);
    set_enabled(false);
  }

  bool open(struct ged* gedp, const std::string& object_name)
  {
    if (dirty_) {
      return fail("Save or discard the current combination changes before opening another object");
    }
    ged_ = gedp;
    object_name_ = object_name;
    if (!load()) {
      ged_ = nullptr;
      object_name_.clear();
      set_enabled(false);
      return false;
    }
    set_enabled(true);
    dirty_ = false;
    update_buttons();
    owner_.copy_label(("Combination Editor - " + object_name_).c_str());
    return true;
  }

  bool save()
  {
    if (!ged_ || !ged_->dbip || !directory_) {
      return fail("No combination is open for saving");
    }
    if (!dirty_) {
      return true;
    }
    if (ged_->dbip->dbi_read_only || db_version(ged_->dbip) < 5) {
      return fail("The database is read-only or does not support editing");
    }
    struct directory* current =
      db_lookup(ged_->dbip, object_name_.c_str(), LOOKUP_QUIET);
    if (current != directory_ || current->d_addr != saved_address_ ||
        current->d_len != saved_length_ ||
        current->d_minor_type != ID_COMBINATION) {
      return fail("The combination changed outside this editor; reload it before saving");
    }
    const char* tree_value = tree_->value() ? tree_->value() : "";
    const char* shader_value = shader_->value() ? shader_->value() : "";
    const char* material_value = material_->value() ? material_->value() : "";
    if (std::strlen(tree_value) > kMaxTreeBytes) {
      return fail("The Boolean tree exceeds the editor size limit");
    }
    if (std::strlen(shader_value) > kMaxTextBytes ||
        std::strlen(material_value) > kMaxTextBytes) {
      return fail("Shader or material text exceeds the editor size limit");
    }
    if (!valid_color(color_->value())) {
      return fail("Color must be empty or three whole numbers from 0 through 255");
    }

    long region_id = 0;
    long air = 0;
    long los = 0;
    long material_id = 0;
    if (region_->value() &&
        (!parse_integer(region_id_->value(), region_id) ||
         !parse_integer(air_->value(), air) ||
         !parse_integer(los_->value(), los) ||
         !parse_integer(material_id_->value(), material_id))) {
      return fail("Region identifiers must be whole numbers");
    }

    struct rt_db_internal replacement;
    RT_DB_INTERNAL_INIT(&replacement);
    if (rt_db_get_internal(&replacement, current, ged_->dbip, nullptr) < 0) {
      return fail("Unable to reload the combination for validation");
    }

    std::vector<std::string> storage;
    storage.reserve(20);
    const auto add = [&storage](const char* key, std::string value) {
      storage.emplace_back(key);
      storage.push_back(std::move(value));
    };
    add("region", region_->value() ? "yes" : "no");
    if (region_->value()) {
      add("id", integer_text(region_id));
      add("air", integer_text(air));
      add("los", integer_text(los));
      add("GIFTmater", integer_text(material_id));
    }
    add("rgb", color_->value() && *color_->value() ? color_->value() : "none");
    add("shader", *shader_value ? shader_value : "none");
    add("material", *material_value ? material_value : "none");
    add("inherit", inherit_->value() ? "yes" : "no");
    add("tree", *tree_value ? tree_value : "none");

    std::vector<const char*> arguments;
    arguments.reserve(storage.size());
    for (const auto& value : storage) {
      arguments.push_back(value.c_str());
    }
    struct bu_vls errors = BU_VLS_INIT_ZERO;
    const int adjusted = OBJ[ID_COMBINATION].ft_adjust(
      &errors, &replacement, static_cast<int>(arguments.size()),
      arguments.data());
    if (adjusted & BRLCAD_ERROR) {
      const std::size_t length = bu_vls_strlen(&errors);
      const std::size_t accepted = std::min(length, kMaxTextBytes);
      std::string message = length
                              ? std::string(bu_vls_cstr(&errors), accepted)
                              : "Invalid combination values";
      if (accepted != length) {
        message += "\n[validation diagnostic truncated]";
      }
      bu_vls_free(&errors);
      rt_db_free_internal(&replacement);
      return fail(message);
    }
    bu_vls_free(&errors);
    if (rt_db_put_internal(current, ged_->dbip, &replacement) < 0) {
      rt_db_free_internal(&replacement);
      return fail("BRL-CAD could not write the combination");
    }
    saved_address_ = current->d_addr;
    saved_length_ = current->d_len;
    dirty_ = false;
    update_buttons();
    if (changed_handler_) {
      changed_handler_();
    }
    last_error_.clear();
    return true;
  }

  void reload()
  {
    if (dirty_ &&
        fl_choice("Discard unsaved changes to %s?", "Cancel", "Discard",
                  nullptr, object_name_.c_str()) != 1) {
      return;
    }
    if (!load()) {
      fl_alert("%s", last_error_.c_str());
    }
  }

  void resize(int width, int height)
  {
    heading_->resize(kMargin, kMargin, width - 2 * kMargin, kRowHeight);
    const int input_width = std::max(1, width - kLabelWidth - 2 * kMargin);
    Fl_Input* inputs[] = {region_id_, air_, los_, material_id_,
                          color_, shader_, material_};
    for (int index = 0; index < 7; ++index) {
      inputs[index]->resize(kMargin + kLabelWidth,
                            kMargin + (index + 2) * kRowHeight,
                            input_width, kRowHeight);
    }
    tree_->resize(kMargin + kLabelWidth, kMargin + 9 * kRowHeight,
                  input_width,
                  std::max(1, height - 10 * kRowHeight -
                                kButtonHeight - 3 * kMargin));
    save_->position(width - 3 * (kButtonWidth + kColumnGap),
                    height - kButtonHeight - kMargin);
    reload_->position(width - 2 * (kButtonWidth + kColumnGap),
                      height - kButtonHeight - kMargin);
    close_->position(width - kButtonWidth - kMargin,
                     height - kButtonHeight - kMargin);
    add_union_->position(kMargin, height - kButtonHeight - kMargin);
    add_subtract_->position(kMargin + kButtonWidth + kColumnGap,
                            height - kButtonHeight - kMargin);
    add_intersect_->position(kMargin + 2 * (kButtonWidth + kColumnGap),
                             height - kButtonHeight - kMargin);
  }

  const std::string& object_name() const noexcept { return object_name_; }
  const std::string& last_error() const noexcept { return last_error_; }
  bool dirty() const noexcept { return dirty_; }
  void changed_handler(ChangedHandler handler)
  {
    changed_handler_ = std::move(handler);
  }
  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

private:
  Fl_Input* make_input(const char* label, int row)
  {
    auto* input = new Fl_Input(kMargin + kLabelWidth,
                               kMargin + row * kRowHeight,
                               owner_.w() - kLabelWidth - 2 * kMargin,
                               kRowHeight, label);
    input->align(FL_ALIGN_LEFT);
    return input;
  }

  bool load()
  {
    if (!ged_ || !ged_->dbip || object_name_.empty()) {
      return fail("Open a database and select a combination first");
    }
    struct directory* current =
      db_lookup(ged_->dbip, object_name_.c_str(), LOOKUP_QUIET);
    if (!current || (current->d_flags & RT_DIR_COMB) == 0 ||
        current->d_minor_type != ID_COMBINATION) {
      return fail("The selected object is not a combination");
    }
    struct rt_db_internal internal;
    RT_DB_INTERNAL_INIT(&internal);
    if (rt_db_get_internal(&internal, current, ged_->dbip, nullptr) < 0) {
      return fail("Unable to read the selected combination");
    }
    auto* combination = static_cast<struct rt_comb_internal*>(internal.idb_ptr);
    RT_CK_COMB(combination);
    struct bu_vls tree = BU_VLS_INIT_ZERO;
    db_tree_list(&tree, combination->tree);
    if (bu_vls_strlen(&tree) > kMaxTreeBytes) {
      bu_vls_free(&tree);
      rt_db_free_internal(&internal);
      return fail("The combination tree exceeds the editor size limit");
    }
    if (bu_vls_strlen(&combination->shader) > kMaxTextBytes ||
        bu_vls_strlen(&combination->material) > kMaxTextBytes) {
      bu_vls_free(&tree);
      rt_db_free_internal(&internal);
      return fail("The combination shader or material exceeds the editor size limit");
    }

    loading_ = true;
    directory_ = current;
    saved_address_ = current->d_addr;
    saved_length_ = current->d_len;
    region_->value(combination->region_flag != 0);
    inherit_->value(combination->inherit != 0);
    region_id_->value(integer_text(combination->region_id).c_str());
    air_->value(integer_text(combination->aircode).c_str());
    los_->value(integer_text(combination->los).c_str());
    material_id_->value(integer_text(combination->GIFTmater).c_str());
    if (combination->rgb_valid) {
      color_->value((integer_text(combination->rgb[0]) + " " +
                     integer_text(combination->rgb[1]) + " " +
                     integer_text(combination->rgb[2]))
                      .c_str());
    } else {
      color_->value("");
    }
    shader_->value(bu_vls_cstr(&combination->shader));
    material_->value(bu_vls_cstr(&combination->material));
    tree_->value(bu_vls_cstr(&tree));
    heading_->copy_label(("Combination: " + object_name_).c_str());
    bu_vls_free(&tree);
    rt_db_free_internal(&internal);
    loading_ = false;
    dirty_ = false;
    update_buttons();
    last_error_.clear();
    return true;
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* widgets[] = {region_, region_id_, air_, los_, material_id_,
                            color_, shader_, material_, inherit_, tree_, save_,
                            reload_, add_union_, add_subtract_, add_intersect_};
    for (auto* widget : widgets) {
      if (enabled) {
        widget->activate();
      } else {
        widget->deactivate();
      }
    }
    update_buttons();
  }

  void update_buttons()
  {
    if (ged_ && dirty_) {
      save_->activate();
    } else {
      save_->deactivate();
    }
  }

  void mark_modified()
  {
    if (!loading_ && ged_) {
      dirty_ = true;
      update_buttons();
    }
  }

  void add_member(Fl_Widget* source)
  {
    if (!ged_ || !ged_->dbip) {
      fail("Open a combination before adding a member");
      return;
    }
    const char* entered = fl_input("Database object to add as a member:", "");
    if (!entered) {
      return;
    }
    const std::string member(entered);
    if (member.empty() || member == object_name_ ||
        !db_lookup(ged_->dbip, member.c_str(), LOOKUP_QUIET)) {
      fail("Enter the name of an existing object other than this combination");
      return;
    }

    const char* operation = source == add_subtract_
                              ? "-"
                              : (source == add_intersect_ ? "n" : "u");
    struct bu_vls encoded = BU_VLS_INIT_ZERO;
    bu_vls_encode(&encoded, member.c_str());
    const std::string leaf = "l " + std::string(bu_vls_cstr(&encoded));
    bu_vls_free(&encoded);
    const std::string current = tree_->value() ? tree_->value() : "";
    const std::string replacement =
      current.empty() || current == "none"
        ? leaf
        : std::string(operation) + " {" + current + "} {" + leaf + "}";
    if (replacement.size() > kMaxTreeBytes) {
      fail("Adding that member would exceed the Boolean tree size limit");
      return;
    }
    tree_->value(replacement.c_str());
    mark_modified();
  }

  bool request_close()
  {
    if (!dirty_) {
      owner_.hide();
      return true;
    }
    const int decision = fl_choice(
      "Save changes to %s before closing?", "Cancel", "Discard", "Save",
      object_name_.c_str());
    if (decision == 0) {
      return false;
    }
    if (decision == 2 && !save()) {
      fl_alert("%s", last_error_.c_str());
      return false;
    }
    if (decision == 1) {
      dirty_ = false;
      update_buttons();
    }
    owner_.hide();
    return true;
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return false;
  }

  static void save_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.save()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }
  static void reload_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->reload();
  }
  static void close_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->request_close();
  }
  static void modified_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->mark_modified();
  }
  static void add_member_callback(Fl_Widget* widget, void* data)
  {
    static_cast<Impl*>(data)->add_member(widget);
  }

  Fl_Ged_Combination_Editor& owner_;
  Fl_Box* heading_ = nullptr;
  Fl_Check_Button* region_ = nullptr;
  Fl_Check_Button* inherit_ = nullptr;
  Fl_Input* region_id_ = nullptr;
  Fl_Input* air_ = nullptr;
  Fl_Input* los_ = nullptr;
  Fl_Input* material_id_ = nullptr;
  Fl_Input* color_ = nullptr;
  Fl_Input* shader_ = nullptr;
  Fl_Input* material_ = nullptr;
  Fl_Multiline_Input* tree_ = nullptr;
  Fl_Button* add_union_ = nullptr;
  Fl_Button* add_subtract_ = nullptr;
  Fl_Button* add_intersect_ = nullptr;
  Fl_Button* save_ = nullptr;
  Fl_Button* reload_ = nullptr;
  Fl_Button* close_ = nullptr;
  struct ged* ged_ = nullptr;
  struct directory* directory_ = nullptr;
  b_off_t saved_address_ = 0;
  size_t saved_length_ = 0;
  std::string object_name_;
  std::string last_error_;
  bool dirty_ = false;
  bool loading_ = false;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_Combination_Editor::Fl_Ged_Combination_Editor()
  : Fl_Double_Window(kInitialWidth, kInitialHeight, "Combination Editor"),
    impl_(std::make_unique<Impl>(*this))
{
}

Fl_Ged_Combination_Editor::~Fl_Ged_Combination_Editor() = default;
bool Fl_Ged_Combination_Editor::open(struct ged* gedp,
                                     const std::string& object_name)
{
  return impl_->open(gedp, object_name);
}
bool Fl_Ged_Combination_Editor::save() { return impl_->save(); }
void Fl_Ged_Combination_Editor::reload() { impl_->reload(); }
bool Fl_Ged_Combination_Editor::dirty() const noexcept
{
  return impl_->dirty();
}
const std::string& Fl_Ged_Combination_Editor::object_name() const noexcept
{
  return impl_->object_name();
}
const std::string& Fl_Ged_Combination_Editor::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Combination_Editor::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler(std::move(handler));
}
void Fl_Ged_Combination_Editor::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
}
void Fl_Ged_Combination_Editor::resize(int x, int y, int width, int height)
{
  Fl_Double_Window::resize(x, y, width, height);
  if (impl_) {
    impl_->resize(width, height);
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
