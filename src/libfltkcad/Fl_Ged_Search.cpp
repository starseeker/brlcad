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

#include "fltkcad/Fl_Ged_Search.H"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Input.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <cctype>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "ged.h"
#include "raytrace.h"
}

namespace {

constexpr int kWindowWidth = 620;
constexpr int kWindowHeight = 520;
constexpr int kMargin = 10;
constexpr int kRowHeight = 28;
constexpr int kGap = 6;
constexpr int kTypeWidth = 145;
constexpr int kButtonWidth = 80;
constexpr int kQueryLimit = 1024;
constexpr std::size_t kResultLimit = 20000;
constexpr std::size_t kInitialResultCapacity = 1024;

std::string lower_copy(const char* value)
{
  std::string result = value ? value : "";
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char character) {
                   return static_cast<char>(std::tolower(character));
                 });
  return result;
}

std::string quoted_argument(const std::string& value)
{
  std::string result("\"");
  for (char character : value) {
    if (character == '\\' || character == '"') {
      result.push_back('\\');
    }
    result.push_back(character);
  }
  result.push_back('"');
  return result;
}

} // namespace

class Fl_Ged_Search::Impl
{
public:
  explicit Impl(Fl_Ged_Search& owner) : owner_(owner)
  {
    owner_.begin();
    query_ = new Fl_Input(kMargin + 54, kMargin,
                          owner_.w() - kTypeWidth - 4 * kMargin - 54,
                          kRowHeight, "Name");
    query_->maximum_size(kQueryLimit);
    type_ = new Fl_Choice(owner_.w() - kTypeWidth - kMargin, kMargin,
                          kTypeWidth, kRowHeight, "");
    type_->add("All objects|Primitives|Combinations|Regions");
    type_->value(0);
    search_ = new Fl_Button(kMargin, kMargin + kRowHeight + kGap,
                            kButtonWidth, kRowHeight, "Search");
    count_ = new Fl_Box(kMargin + kButtonWidth + kGap,
                        kMargin + kRowHeight + kGap,
                        owner_.w() - kButtonWidth - 3 * kMargin,
                        kRowHeight, "No database open");
    count_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    results_ = new Fl_Hold_Browser(
      kMargin, kMargin + 2 * (kRowHeight + kGap), owner_.w() - 2 * kMargin,
      owner_.h() - 3 * kRowHeight - 4 * kMargin - 2 * kGap);
    draw_ = new Fl_Button(owner_.w() - 3 * (kButtonWidth + kGap) - kMargin,
                          owner_.h() - kRowHeight - kMargin, kButtonWidth,
                          kRowHeight, "Draw");
    erase_ = new Fl_Button(owner_.w() - 2 * (kButtonWidth + kGap) - kMargin,
                           owner_.h() - kRowHeight - kMargin, kButtonWidth,
                           kRowHeight, "Erase");
    close_ = new Fl_Button(owner_.w() - kButtonWidth - kMargin,
                           owner_.h() - kRowHeight - kMargin, kButtonWidth,
                           kRowHeight, "Close");
    search_->callback(search_callback, this);
    query_->when(FL_WHEN_ENTER_KEY_ALWAYS);
    query_->callback(search_callback, this);
    type_->callback(search_callback, this);
    results_->callback(selection_callback, this);
    draw_->callback(command_callback, this);
    erase_->callback(command_callback, this);
    close_->callback(close_callback, this);
    owner_.end();
    owner_.resizable(results_);
    set_enabled(false);
  }

  void ged(struct ged* value)
  {
    ged_ = value;
    set_enabled(ged_ && ged_->dbip);
    search();
  }

  struct ged* ged() const noexcept { return ged_; }

  bool search()
  {
    results_->clear();
    names_.clear();
    selected_.clear();
    update_action_buttons();
    if (!ged_ || !ged_->dbip) {
      count_->copy_label("No database open");
      return false;
    }
    const std::string query = lower_copy(query_->value());
    names_.reserve(kInitialResultCapacity);
    bool truncated = false;
    struct directory* entry = nullptr;
    FOR_ALL_DIRECTORY_START(entry, ged_->dbip)
      if (!entry->d_namep || entry->d_addr == RT_DIR_PHONY_ADDR ||
          (entry->d_flags & RT_DIR_HIDDEN) || !matches_type(*entry)) {
        continue;
      }
      const std::string name(entry->d_namep);
      if (!query.empty() && lower_copy(name.c_str()).find(query) ==
                              std::string::npos) {
        continue;
      }
      if (names_.size() == kResultLimit) {
        truncated = true;
        continue;
      }
      names_.push_back(name);
    FOR_ALL_DIRECTORY_END;
    std::sort(names_.begin(), names_.end());
    for (const auto& name : names_) {
      results_->add(name.c_str());
    }
    std::string summary = std::to_string(names_.size()) + " result";
    if (names_.size() != 1) {
      summary += 's';
    }
    if (truncated) {
      summary += " (limit reached)";
      fail("Search results exceeded the configured limit");
    } else {
      last_error_.clear();
    }
    count_->copy_label(summary.c_str());
    return !truncated;
  }

  const std::string& selected_object() const noexcept { return selected_; }
  const std::string& last_error() const noexcept { return last_error_; }
  void selection_handler(SelectionHandler handler)
  {
    selection_handler_ = std::move(handler);
  }
  void command_handler(CommandHandler handler)
  {
    command_handler_ = std::move(handler);
  }
  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

  void resize(int width, int height)
  {
    query_->resize(kMargin + 54, kMargin,
                   std::max(1, width - kTypeWidth - 4 * kMargin - 54),
                   kRowHeight);
    type_->position(width - kTypeWidth - kMargin, kMargin);
    count_->resize(kMargin + kButtonWidth + kGap,
                   kMargin + kRowHeight + kGap,
                   std::max(1, width - kButtonWidth - 3 * kMargin),
                   kRowHeight);
    results_->resize(kMargin, kMargin + 2 * (kRowHeight + kGap),
                     std::max(1, width - 2 * kMargin),
                     std::max(1, height - 3 * kRowHeight - 4 * kMargin -
                                   2 * kGap));
    draw_->position(width - 3 * (kButtonWidth + kGap) - kMargin,
                    height - kRowHeight - kMargin);
    erase_->position(width - 2 * (kButtonWidth + kGap) - kMargin,
                     height - kRowHeight - kMargin);
    close_->position(width - kButtonWidth - kMargin,
                     height - kRowHeight - kMargin);
  }

private:
  bool matches_type(const struct directory& entry) const
  {
    switch (type_->value()) {
      case 1: return (entry.d_flags & RT_DIR_SOLID) != 0;
      case 2: return (entry.d_flags & RT_DIR_COMB) != 0;
      case 3: return (entry.d_flags & RT_DIR_REGION) != 0;
      default: return true;
    }
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* widgets[] = {query_, type_, search_, results_};
    for (auto* widget : widgets) {
      enabled ? widget->activate() : widget->deactivate();
    }
    update_action_buttons();
  }

  void update_selection()
  {
    const int index = results_->value() - 1;
    if (index < 0 || static_cast<std::size_t>(index) >= names_.size()) {
      selected_.clear();
    } else {
      selected_ = names_[static_cast<std::size_t>(index)];
      if (selection_handler_) {
        selection_handler_(selected_);
      }
    }
    update_action_buttons();
  }

  void update_action_buttons()
  {
    const bool enabled = ged_ && ged_->dbip && !selected_.empty();
    enabled ? draw_->activate() : draw_->deactivate();
    enabled ? erase_->activate() : erase_->deactivate();
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return false;
  }

  static void search_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->search();
  }
  static void selection_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->update_selection();
  }
  static void command_callback(Fl_Widget* widget, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.command_handler_ || self.selected_.empty()) {
      return;
    }
    const char* command = widget == self.draw_ ? "draw " : "erase ";
    self.command_handler_(command + quoted_argument(self.selected_));
  }
  static void close_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->owner_.hide();
  }

  Fl_Ged_Search& owner_;
  Fl_Input* query_ = nullptr;
  Fl_Choice* type_ = nullptr;
  Fl_Button* search_ = nullptr;
  Fl_Box* count_ = nullptr;
  Fl_Hold_Browser* results_ = nullptr;
  Fl_Button* draw_ = nullptr;
  Fl_Button* erase_ = nullptr;
  Fl_Button* close_ = nullptr;
  struct ged* ged_ = nullptr;
  std::vector<std::string> names_;
  std::string selected_;
  std::string last_error_;
  SelectionHandler selection_handler_;
  CommandHandler command_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_Search::Fl_Ged_Search()
  : Fl_Double_Window(kWindowWidth, kWindowHeight, "Database Search"),
    impl_(std::make_unique<Impl>(*this))
{
}
Fl_Ged_Search::~Fl_Ged_Search() = default;
void Fl_Ged_Search::ged(struct ged* value) { impl_->ged(value); }
struct ged* Fl_Ged_Search::ged() const noexcept { return impl_->ged(); }
bool Fl_Ged_Search::search() { return impl_->search(); }
const std::string& Fl_Ged_Search::selected_object() const noexcept
{
  return impl_->selected_object();
}
const std::string& Fl_Ged_Search::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Search::selection_handler(SelectionHandler handler)
{
  impl_->selection_handler(std::move(handler));
}
void Fl_Ged_Search::command_handler(CommandHandler handler)
{
  impl_->command_handler(std::move(handler));
}
void Fl_Ged_Search::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
}
void Fl_Ged_Search::resize(int x, int y, int width, int height)
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
