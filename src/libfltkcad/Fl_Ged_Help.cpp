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

#include "fltkcad/Fl_Ged_Help.H"

#include "fltkcad/Fl_Litehtml_View.H"
#include "Uri.h"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>
#include <cctype>
#include <cstdlib>

extern "C" {
#include "bu/app.h"
#include "bu/file.h"
}

namespace {

constexpr int kToolbarHeight = 32;
constexpr int kToolbarGap = 4;
constexpr int kSearchWidth = 220;
constexpr int kIndexWidth = 220;
constexpr int kNavigationWidth = 48;
constexpr std::size_t kMaxNavigationHistory = 256;
constexpr const char* kManualSections = "135n";
constexpr int kButtonWidth = 64;
constexpr std::size_t kMaxNavigationPathBytes = 64U * 1024U;

std::string default_documentation()
{
  const char* root = bu_dir(nullptr, 0, BU_DIR_DOC, "html", "manuals", "mged",
                            nullptr);
  if (!root || !*root) {
    return {};
  }
  const std::array<const char*, 2> candidates{{"ged.html", "contents.html"}};
  for (const char* candidate : candidates) {
    const std::filesystem::path path = std::filesystem::path(root) / candidate;
    if (bu_file_exists(path.string().c_str(), nullptr)) {
      return path.string();
    }
  }
  return {};
}

} // namespace

class Fl_Ged_Help::Impl
{
public:
  explicit Impl(Fl_Ged_Help& owner)
    : owner_(owner)
  {
    owner_.begin();
    back_ = new Fl_Button(kToolbarGap, kToolbarGap, kNavigationWidth,
                           kToolbarHeight - 2 * kToolbarGap, "Back");
    forward_ = new Fl_Button(kNavigationWidth + 2 * kToolbarGap, kToolbarGap,
                              kNavigationWidth, kToolbarHeight - 2 * kToolbarGap,
                              "Forward");
    home_ = new Fl_Button(2 * kNavigationWidth + 3 * kToolbarGap, kToolbarGap,
                           kButtonWidth, kToolbarHeight - 2 * kToolbarGap, "Manual");
    back_->callback(navigation_callback, this);
    forward_->callback(navigation_callback, this);
    home_->callback(navigation_callback, this);
    const int search_x = 2 * kNavigationWidth + kButtonWidth + 4 * kToolbarGap;
    search_ = new Fl_Input(search_x, kToolbarGap, kSearchWidth,
                           kToolbarHeight - 2 * kToolbarGap);
    search_->when(FL_WHEN_ENTER_KEY_ALWAYS);
    search_->callback(search_callback, this);
    find_ = new Fl_Button(search_x + kSearchWidth + kToolbarGap, kToolbarGap,
                          kButtonWidth, kToolbarHeight - 2 * kToolbarGap, "Find");
    previous_ = new Fl_Button(find_->x() + kButtonWidth + kToolbarGap,
                              kToolbarGap, kButtonWidth,
                              kToolbarHeight - 2 * kToolbarGap, "Prev");
    next_ = new Fl_Button(previous_->x() + kButtonWidth + kToolbarGap,
                          kToolbarGap, kButtonWidth,
                          kToolbarHeight - 2 * kToolbarGap, "Next");
    status_ = new Fl_Box(0, owner_.h() - kToolbarHeight, owner_.w(), kToolbarHeight);
    status_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    find_->callback(search_callback, this);
    previous_->callback(search_callback, this);
    next_->callback(search_callback, this);

    section_ = new Fl_Choice(kToolbarGap, kToolbarHeight + kToolbarGap,
                              kIndexWidth - 2 * kToolbarGap,
                              kToolbarHeight - 2 * kToolbarGap);
    section_->add("All sections|1 - Programs|3 - Libraries|5 - Formats|n - GED commands");
    section_->value(0);
    section_->callback(index_filter_callback, this);
    filter_ = new Fl_Input(kToolbarGap, 2 * kToolbarHeight + kToolbarGap,
                            kIndexWidth - 2 * kToolbarGap,
                            kToolbarHeight - 2 * kToolbarGap);
    filter_->tooltip("Filter manual page names");
    filter_->when(FL_WHEN_CHANGED);
    filter_->callback(index_filter_callback, this);
    index_ = new Fl_Hold_Browser(0, 3 * kToolbarHeight, kIndexWidth,
                                 owner_.h() - 4 * kToolbarHeight);
    index_->callback(index_callback, this);
    view_ = new Fl_Litehtml_View(kIndexWidth, kToolbarHeight,
                                 owner_.w() - kIndexWidth,
                                 owner_.h() - 2 * kToolbarHeight);
    auto limits = view_->limits();
    limits.max_document_bytes = 8U * 1024U * 1024U;
    limits.max_resource_count = 256;
    limits.max_cached_images = 96;
    view_->limits(limits);
    view_->error_handler([this](const std::string& message) {
      last_error_ = message;
      status_->copy_label(message.c_str());
    });
    view_->link_handler([this](const std::string& uri) { navigate(uri); });
    owner_.end();
    owner_.resizable(view_);
    owner_.size_range(720, 400);
    update_navigation();
  }

  bool load(const std::string& path, bool record = true)
  {
    if (path.empty() || path.size() > kMaxNavigationPathBytes ||
        path.find('\0') != std::string::npos) {
      return fail("Documentation path is empty or exceeds the configured limit");
    }
    if (!view_->load(path)) {
      return fail(view_->last_error().empty() ? "Unable to load documentation"
                                              : view_->last_error());
    }
    current_path_ = path;
    if (record) {
      if (!history_.empty()) history_.resize(history_position_ + 1);
      if (history_.empty() || history_.back() != path) {
        history_.push_back(path);
        if (history_.size() > kMaxNavigationHistory) history_.erase(history_.begin());
      }
      history_position_ = history_.size() - 1;
    }
    update_navigation();
    last_error_.clear();
    owner_.copy_label(view_->document_title().empty()
                        ? "MGED Documentation"
                        : view_->document_title().c_str());
    status_->copy_label("");
    return true;
  }

  void update_navigation()
  {
    if (history_position_ > 0) back_->activate(); else back_->deactivate();
    if (history_position_ + 1 < history_.size()) forward_->activate();
    else forward_->deactivate();
  }

  bool travel(bool forward)
  {
    if (history_.empty() || (forward ? history_position_ + 1 >= history_.size()
                                      : history_position_ == 0)) return false;
    const std::size_t position = forward ? history_position_ + 1 : history_position_ - 1;
    if (!load(history_[position], false)) return false;
    history_position_ = position;
    update_navigation();
    return true;
  }

  struct ManualPage { std::string name; std::string path; char section; };

  bool manuals(const std::string& language, char section)
  {
    const bool valid_language = language.size() == 2 &&
      std::all_of(language.begin(), language.end(),
                  [](unsigned char c) { return std::isalpha(c); });
    const std::string locale = valid_language ? language : "en";
    if (section && std::string(kManualSections).find(section) == std::string::npos)
      return fail("Manual section must be 1, 3, 5, or n");
    const std::filesystem::path root = bu_dir(nullptr, 0, BU_DIR_DOC, "html", nullptr);
    pages_.clear();
    for (const char* current = kManualSections; *current; ++current) {
      const std::string directory = std::string("man") + *current;
      std::error_code error;
      std::filesystem::path path = root / locale / directory;
      if (!std::filesystem::is_directory(path, error)) path = root / directory;
      error.clear();
      std::filesystem::directory_iterator entry(path, error), end;
      while (!error && entry != end) {
        if (entry->is_regular_file(error) && entry->path().extension() == ".html")
          pages_.push_back({entry->path().stem().string(), entry->path().string(), *current});
        entry.increment(error);
      }
    }
    std::sort(pages_.begin(), pages_.end(), [](const ManualPage& a, const ManualPage& b) {
      return a.name == b.name ? a.section < b.section : a.name < b.name;
    });
    section_->value(section ? static_cast<int>(std::string(kManualSections).find(section)) + 1 : 0);
    filter_index();
    if (pages_.empty()) return fail("BRL-CAD HTML manual pages are not installed");
    return true;
  }

  void filter_index()
  {
    std::string filter = filter_->value();
    std::transform(filter.begin(), filter.end(), filter.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    index_->clear();
    indexed_pages_.clear();
    for (std::size_t page = 0; page < pages_.size(); ++page) {
      const ManualPage& manual = pages_[page];
      if (section_->value() > 0 &&
          manual.section != kManualSections[section_->value() - 1]) continue;
      std::string name = manual.name;
      std::transform(name.begin(), name.end(), name.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (name.find(filter) == std::string::npos) continue;
      index_->add((manual.name + " (" + manual.section + ")").c_str());
      indexed_pages_.push_back(page);
    }
  }

  static void index_filter_callback(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->filter_index(); }

  static void index_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const int selected = self.index_->value();
    if (selected > 0 && static_cast<std::size_t>(selected) <= self.indexed_pages_.size())
      self.load(self.pages_[self.indexed_pages_[static_cast<std::size_t>(selected - 1)]].path);
  }

  static void navigation_callback(Fl_Widget* widget, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (widget == self.back_) self.travel(false);
    else if (widget == self.forward_) self.travel(true);
    else self.owner_.load_default();
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    status_->copy_label(last_error_.c_str());
    return false;
  }

  void navigate(const std::string& uri)
  {
    if (uri.empty() || uri.size() > kMaxNavigationPathBytes) {
      fail("Documentation link exceeds the configured limit");
      return;
    }
    if (uri.front() == '#') {
      if (!view_->scroll_to_fragment(uri.substr(1))) {
        fail("Documentation section was not found");
      }
      return;
    }
    const auto path = fltkcad::detail::local_path(uri);
    if (!path) {
      fail("External documentation links are disabled");
      return;
    }
    const bool current_document = !current_path_.empty() &&
      fltkcad::detail::file_uri(*path) == fltkcad::detail::file_uri(current_path_);
    if (current_document || load(path->string())) {
      const auto fragment = uri.find('#');
      if (fragment != std::string::npos &&
          !view_->scroll_to_fragment(uri.substr(fragment + 1))) {
        fail("Documentation section was not found");
      }
    }
  }

  void search(Fl_Widget* source)
  {
    if (source == previous_) {
      view_->find_previous();
    } else if (source == next_) {
      view_->find_next();
    } else {
      const char* query = search_->value();
      view_->find(query ? query : "");
    }
    if (view_->match_count() == 0) {
      status_->copy_label("No matches");
      return;
    }
    const std::string message = std::to_string(view_->current_match()) + " of " +
                                std::to_string(view_->match_count());
    status_->copy_label(message.c_str());
  }

  static void search_callback(Fl_Widget* widget, void* data)
  {
    static_cast<Impl*>(data)->search(widget);
  }

  Fl_Ged_Help& owner_;
  Fl_Button* back_ = nullptr;
  Fl_Button* forward_ = nullptr;
  Fl_Button* home_ = nullptr;
  Fl_Choice* section_ = nullptr;
  Fl_Input* filter_ = nullptr;
  Fl_Hold_Browser* index_ = nullptr;
  std::vector<ManualPage> pages_;
  std::vector<std::size_t> indexed_pages_;
  std::vector<std::string> history_;
  std::size_t history_position_ = 0;
  Fl_Input* search_ = nullptr;
  Fl_Button* find_ = nullptr;
  Fl_Button* previous_ = nullptr;
  Fl_Button* next_ = nullptr;
  Fl_Box* status_ = nullptr;
  Fl_Litehtml_View* view_ = nullptr;
  std::string current_path_;
  std::string last_error_;
};

Fl_Ged_Help::Fl_Ged_Help(int width, int height, const char* title)
  : Fl_Double_Window(width, height, title),
    impl_(std::make_unique<Impl>(*this))
{
}

Fl_Ged_Help::~Fl_Ged_Help() = default;

bool Fl_Ged_Help::load_default()
{
  const char* locale = std::getenv("LANG");
  const std::string language = locale && std::char_traits<char>::length(locale) >= 2
                                  ? std::string(locale, 2) : "en";
  impl_->manuals(language, '\0');
  const std::string path = default_documentation();
  if (path.empty()) {
    const auto page = std::find_if(impl_->pages_.begin(), impl_->pages_.end(),
      [](const Impl::ManualPage& manual) { return manual.name == "mged"; });
    if (page != impl_->pages_.end()) return impl_->load(page->path);
    return impl_->fail("MGED HTML documentation is not installed");
  }
  return impl_->load(path);
}

bool Fl_Ged_Help::manuals(const std::string& language, char section)
{ return impl_->manuals(language, section); }
std::size_t Fl_Ged_Help::manual_count() const noexcept { return impl_->pages_.size(); }
bool Fl_Ged_Help::back() { return impl_->travel(false); }
bool Fl_Ged_Help::forward() { return impl_->travel(true); }

bool Fl_Ged_Help::load(const std::string& path) { return impl_->load(path); }

const std::string& Fl_Ged_Help::last_error() const noexcept
{
  return impl_->last_error_;
}

Fl_Ged_About::Fl_Ged_About(const char* backend_name)
  : Fl_Double_Window(460, 250, "About MGED")
{
  begin();
  auto* title = new Fl_Box(20, 18, 420, 34, "MGED");
  title->labelfont(FL_BOLD);
  title->labelsize(24);
  const std::string version =
    std::string("BRL-CAD ") + FLTKCAD_BRLCAD_VERSION;
  auto* version_box = new Fl_Box(20, 54, 420, 24);
  version_box->copy_label(version.c_str());
  version_box->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE);
  auto* description = new Fl_Box(
    30, 84, 400, 70,
    "Multi-device Geometry Editor\n"
    "An interactive solid-modeling application in the BRL-CAD package.");
  description->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE | FL_ALIGN_WRAP);
  const std::string backend =
    std::string("Display backend: ") +
    (backend_name && *backend_name ? backend_name : "not initialized");
  auto* details = new Fl_Box(20, 158, 420, 24);
  details->copy_label(backend.c_str());
  details->align(FL_ALIGN_CENTER | FL_ALIGN_INSIDE);
  auto* docs = new Fl_Button(126, 202, 96, 30, "Documentation");
  auto* close = new Fl_Button(238, 202, 96, 30, "Close");
  docs->callback(docs_callback, this);
  close->callback(close_callback, this);
  end();
  set_modal();
}

Fl_Ged_About::~Fl_Ged_About() = default;

void Fl_Ged_About::close_callback(Fl_Widget*, void* data)
{
  static_cast<Fl_Ged_About*>(data)->hide();
}

void Fl_Ged_About::docs_callback(Fl_Widget*, void* data)
{
  auto& self = *static_cast<Fl_Ged_About*>(data);
  if (!self.help_) {
    self.help_ = std::make_unique<Fl_Ged_Help>();
  }
  if (!self.help_->load_default()) {
    fl_alert("%s", self.help_->last_error().c_str());
    return;
  }
  self.hide();
  self.help_->show();
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
