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
#include "bu/str.h"

#include "FltkGedWindow.H"

#include "fltkcad/Fl_Cad_Terminal.H"
#include "fltkcad/Fl_Cad_Quad_View.H"
#include "fltkcad/Fl_Ged_Accessibility.H"
#include "fltkcad/Fl_Ged_Attributes.H"
#include "fltkcad/Fl_Ged_Combination_Editor.H"
#include "fltkcad/Fl_Ged_Create_Primitive.H"
#include "fltkcad/Fl_Ged_Help.H"
#include "fltkcad/Fl_Ged_Hierarchy.H"
#include "fltkcad/Fl_Ged_Interaction.H"
#include "Fl_Ged_Localization.H"
#include "fltkcad/Fl_Ged_Primitive_Editor.H"
#include "fltkcad/Fl_Ged_Preferences.H"
#include "fltkcad/Fl_Ged_Render_Tools.H"
#include "fltkcad/Fl_Ged_Search.H"
#include "fltkcad/Fl_Ged_Sketch_Editor.H"
#include "fltkcad/Fl_Ged_Specialist_Tools.H"
#include "fltkcad/Fl_Ged_Tcl.H"
#include "fltkcad/Fl_Ged_Transform_Editor.H"
#include "fltkcad/Fl_Ged_View_Controls.H"

#include <FL/Fl.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <FL/Fl_Preferences.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Tile.H>
#include <FL/Fl_Text_Display.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iterator>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/process.h"
#include "bu/vls.h"
#include "bv/util.h"
#include "bn/mat.h"
#include "dm/faceplate.h"
#include "dm/fbserv.h"
#include "ged.h"
#include "ged/edit.h"
}

namespace {

constexpr int kMenuHeight = 25;
constexpr int kPanelWidth = 420;
constexpr int kPanelHeight = 600;
enum class FaceplateAction {
  Title, Oblique, Top, Right, Front, Isometric, Restore, Save,
  Adc, Sliders, ZoomIn, ZoomOut
};
struct FaceplateMenuItem {
  const char* label;
  FaceplateAction action;
};
constexpr FaceplateMenuItem kFaceplateMenu[]{
  {"BUTTON MENU", FaceplateAction::Title},
  {"35,25", FaceplateAction::Oblique}, {"Top", FaceplateAction::Top},
  {"Right", FaceplateAction::Right}, {"Front", FaceplateAction::Front},
  {"45,45", FaceplateAction::Isometric}, {"Restore View", FaceplateAction::Restore},
  {"Save View", FaceplateAction::Save}, {"Ang/Dist Curs", FaceplateAction::Adc},
  {"Sliders", FaceplateAction::Sliders}, {"Zoom In 2X", FaceplateAction::ZoomIn},
  {"Zoom Out 2X", FaceplateAction::ZoomOut}
};
struct EditMenuAction {
  const char* label;
  const char* command;
};
constexpr EditMenuAction kEditMenuActions[]{
  {"Accept", "accept"}, {"Apply", "sed_apply"},
  {"Reset", "sed_reset"}, {"Reject", "reject"}
};
constexpr const char* kSliderLabels[]{"xslew", "yslew", "zslew", "scale", "xrot", "yrot", "zrot"};
constexpr int kInitialViewPercent = 68;
constexpr std::size_t kMaxPendingLogBytes = 4U * 1024U * 1024U;
constexpr std::size_t kMaxCommandArguments = 4096;
constexpr std::size_t kMaxPendingCommandBytes = 256U * 1024U;
constexpr std::size_t kMaxPromptBytes = 1024;
constexpr const char* kDefaultPrompt = "mged> ";
constexpr const char* kContinuationPrompt = "...> ";
constexpr const char* kOutputTruncation =
  "\n[GED output truncated by application limit]\n";
constexpr std::size_t kMaxStoredViews = 10;
constexpr const char* kPreferencesVendor = "brlcad.org";
constexpr const char* kPreferencesApplication = "fltkged";
constexpr int kMaximumPreferenceHistory = 100000;
constexpr int kMaximumPreferenceScrollback = 1000000;
constexpr double kProcessPollIntervalSeconds = 0.02;

const std::set<std::string> kDatabaseCommands{
  "adjust", "arced", "attr", "bo", "bot", "clone", "comb", "concat", "cp",
  "dbconcat", "dbupgrade", "facetize", "g", "in", "kill", "killall",
  "killtree", "make", "mater", "move", "mv", "orotate", "oscale",
  "otranslate", "put", "r", "rarc", "region", "rm", "shader", "title",
  "units"
};

bool is_space(char value)
{
  return std::isspace(static_cast<unsigned char>(value)) != 0;
}

std::string quote(const std::string& value)
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

std::vector<std::string> parse_command(const std::string& command)
{
  std::vector<char> text(command.begin(), command.end());
  text.push_back('\0');
  std::vector<char*> argv(command.size() + 1, nullptr);
  const int argc = bu_argv_from_string(argv.data(), command.size(), text.data());

  std::vector<std::string> result;
  result.reserve(std::max(0, argc));
  for (int index = 0; index < argc; ++index) {
    result.emplace_back(argv[static_cast<std::size_t>(index)]);
  }
  return result;
}

std::size_t argument_bytes(const std::vector<std::string>& arguments)
{
  std::size_t total = 0;
  for (const auto& argument : arguments) {
    const auto maximum = std::numeric_limits<std::size_t>::max();
    if (total == maximum || argument.size() > maximum - total - 1) {
      return std::numeric_limits<std::size_t>::max();
    }
    total += argument.size() + 1;
  }
  return total;
}

std::string display_path(const std::vector<std::string>& components)
{
  std::string path;
  for (const auto& component : components) {
    path.push_back('/');
    path.append(component);
  }
  return path;
}

bool may_change_database(const std::vector<std::string>& arguments)
{
  return !arguments.empty() &&
         kDatabaseCommands.find(arguments.front()) != kDatabaseCommands.end();
}

const char* backend_name(Fl_Cad_View::Backend backend)
{
  switch (backend) {
    case Fl_Cad_View::Backend::SystemOpenGL: return "System OpenGL";
    case Fl_Cad_View::Backend::OSMesa: return "OSMesa";
    case Fl_Cad_View::Backend::Auto: return "Automatic";
  }
  return "Unknown";
}

} // namespace

class FltkGedWindow::Impl
{
public:
  struct ViewSnapshot
  {
    point_t center{};
    vect_t aet{};
    double size = 0.0;
    double perspective = 0.0;
  };

  struct Faceplate
  {
    Impl* owner = nullptr;
    Fl_Cad_View* widget = nullptr;
    struct dm_faceplate_controls menu{};
    struct dm_faceplate_controls sliders{};
    struct dm_faceplate_gui gui{};
    bool enabled = false;
    bool sliders_enabled = true;
    std::array<double, std::size(kSliderLabels)> values{};
  };

  struct IoListener
  {
    Impl* owner = nullptr;
    struct ged_subprocess* process = nullptr;
    bu_process_io_t type = BU_PROCESS_STDOUT;
    ged_io_func_t callback = nullptr;
    void* callback_data = nullptr;
    int descriptor = -1;
    bool watching = false;
    bool polling = false;
  };

  struct FbWatch
  {
    Impl* owner = nullptr;
    struct fbserv_obj* server = nullptr;
    int client_index = -1;
    int descriptor = -1;
    bool watching = false;
  };

  Impl(FltkGedWindow& owner, Fl_Cad_View::Backend backend)
    : owner_(owner)
  {
    ged_ = ged_create();
    tcl_ = std::make_unique<Fl_Ged_Tcl>();
    if (!tcl_->available()) command_language_ = FltkGedWindow::CommandLanguage::Ged;

    owner_.begin();
    menu_ = new Fl_Menu_Bar(0, 0, owner_.w(), kMenuHeight);
    menu_->add(FLTKGED_TR("&File/&Open..."), FL_CTRL + 'o', open_menu, this);
    menu_->add(FLTKGED_TR("&File/&New Database..."), FL_CTRL + FL_SHIFT + 'n',
               new_database_menu, this);
    menu_->add(FLTKGED_TR("&File/&Import Database..."), 0,
               import_database_menu, this);
    menu_->add(FLTKGED_TR("&File/&Export Selected..."), 0,
               export_selected_menu, this);
    menu_->add(FLTKGED_TR("&File/Load &Tcl Script..."), 0, source_tcl_menu,
               this);
    menu_->add(FLTKGED_TR("&File/Enable Legacy T&k Scripts"), 0,
               enable_tk_menu, this);
    menu_->add(FLTKGED_TR("&File/&Quit"), FL_CTRL + 'q', quit_menu, this);
    menu_->add(FLTKGED_TR("&Edit/&Refresh Database Panels"), FL_CTRL + 'l',
               refresh_menu, this);
    menu_->add(FLTKGED_TR("&Edit/&Search Database..."), FL_CTRL + 'g',
               search_menu, this);
    menu_->add(FLTKGED_TR("&Edit/Edit Selected &Primitive..."), FL_CTRL + 'e',
               edit_primitive_menu, this);
    menu_->add(FLTKGED_TR("&Edit/Edit Selected &Combination..."), 0,
               edit_combination_menu, this);
    menu_->add(FLTKGED_TR("&Edit/Edit Selected &Sketch..."), 0,
               edit_sketch_menu, this);
    menu_->add(FLTKGED_TR("&Edit/Object and Member &Transform..."),
               FL_CTRL + 'm',
               transform_editor_menu, this);
    menu_->add(FLTKGED_TR("&Create/&Primitive..."), FL_CTRL + 'n',
               create_primitive_menu, this);
    menu_->add(FLTKGED_TR("&View/&Single View"), 0, single_view_menu, this);
    menu_->add(FLTKGED_TR("&View/&Quad View"), 0, quad_view_menu, this);
    menu_->add(FLTKGED_TR("&View/&Redraw"), FL_CTRL + 'r', redraw_menu, this);
    menu_->add(FLTKGED_TR("&View/&Autoview"), 0, autoview_menu, this);
    menu_->add(FLTKGED_TR("&View/Top"), 't', top_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Front"), 'f', front_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Right"), 'r', right_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Left"), 'l', left_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Azimuth 35, Elevation 25"), 0, aet_view_menu,
               this);
    menu_->add(FLTKGED_TR("&View/Stored Views/&Add Current"), 0,
               add_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Stored Views/&Next"), FL_CTRL + FL_Right,
               next_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Stored Views/&Previous"), FL_CTRL + FL_Left,
               previous_view_menu, this);
    menu_->add(FLTKGED_TR("&View/Stored Views/&Delete Current"), 0,
               delete_view_menu, this);
    menu_->add(FLTKGED_TR("&Settings/&Preferences..."), 0, preferences_menu, this);
    menu_->add(FLTKGED_TR("&Settings/Command Language/&Tcl"), 0,
               tcl_prompt_menu, this, FL_MENU_RADIO | FL_MENU_VALUE);
    menu_->add(FLTKGED_TR("&Settings/Command Language/&Raw libged"), 0,
               ged_prompt_menu, this, FL_MENU_RADIO);
    menu_->add(FLTKGED_TR("&Settings/&Combined Windows"), 0,
               combined_menu, this, FL_MENU_TOGGLE);
    menu_->add(FLTKGED_TR("&Modes/&Faceplate"), FL_F + 7,
               faceplate_menu, this, FL_MENU_TOGGLE | FL_MENU_VALUE);
    menu_->add(FLTKGED_TR("&Modes/Faceplate &GUI"), FL_F + 8,
               faceplate_gui_menu, this, FL_MENU_TOGGLE);
    menu_->add(FLTKGED_TR("&Modes/Mouse/&Default"), 0, navigate_mode_menu, this);
    menu_->add(FLTKGED_TR("&Modes/Mouse/&Pick Object"), 0, select_mode_menu, this);
    menu_->add(FLTKGED_TR("&Modes/Mouse/&Query Ray"), 0, query_mode_menu, this);
    menu_->add(FLTKGED_TR("&Modes/Mouse/Sweep &Zoom Rectangle"), 0, zoom_mode_menu, this);
    menu_->add(FLTKGED_TR("&Modes/Mouse/Sweep &Select Rectangle"), 0, rectangle_mode_menu, this);
    menu_->add(FLTKGED_TR("&Misc/Query &Center Ray"), 0, center_ray_menu,
               this);
    menu_->add(FLTKGED_TR("&Tools/&Raytrace and Framebuffer..."), 0,
               render_tools_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&Abort Raytrace"), 0, abort_raytrace_menu,
               this);
    menu_->add(FLTKGED_TR("&Tools/&Analyze Selected"), 0,
               analyze_selected_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&GQA Selected"), 0, gqa_selected_menu, this);
    menu_->add(FLTKGED_TR("&Tools/Check &Overlaps for Selected"), 0,
               check_overlaps_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&BOT, NMG, and Material Tools..."), 0,
               specialist_tools_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&Database and View Controls..."), 0,
               panels_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&Command Window"), 0, command_window_menu, this);
    menu_->add(FLTKGED_TR("&Tools/&Graphics Window"), 0, graphics_window_menu, this);
    menu_->add(FLTKGED_TR("&Help/&Documentation"), FL_F + 1,
               documentation_menu, this);
    menu_->add(FLTKGED_TR("&Help/&About"), 0, about_menu, this);

    const int content_height = owner_.h() - kMenuHeight;
    main_tile_ = new Fl_Tile(0, kMenuHeight, owner_.w(), content_height);
    main_tile_->begin();
    terminal_ = new Fl_Cad_Terminal(0, kMenuHeight, owner_.w(), content_height);
    main_tile_->end();
    main_tile_->resizable(terminal_);
    owner_.end();
    owner_.resizable(main_tile_);
    owner_.callback(quit_menu, this);
    owner_.copy_label("MGED Command Window - fltkged");

    graphics_window_ = std::make_unique<Fl_Double_Window>(owner_.w(), owner_.w(),
                                                         "MGED Graphics Window - fltkged");
    graphics_window_->begin();
    graphics_menu_ = new Fl_Menu_Bar(0, 0, owner_.w(), kMenuHeight);
    views_ = new Fl_Cad_Quad_View(0, kMenuHeight, owner_.w(),
                                  graphics_window_->h() - kMenuHeight, backend);
    constexpr std::size_t oblique_view = 1;
    views_->active_index(oblique_view);
    graphics_window_->end();
    graphics_window_->resizable(views_);
    graphics_window_->callback(quit_menu, this);
    sync_menus();

    panels_window_ = std::make_unique<Fl_Double_Window>(kPanelWidth, kPanelHeight,
                                                       "Database and View Controls");
    panels_window_->begin();
    panels_ = new Fl_Tabs(0, 0, kPanelWidth, kPanelHeight);
    panels_->begin();
    hierarchy_group_ =
      new Fl_Group(0, kMenuHeight, kPanelWidth,
                   kPanelHeight - kMenuHeight, "Hierarchy");
    hierarchy_group_->begin();
    hierarchy_ = new Fl_Ged_Hierarchy(
      4, kMenuHeight + 4, kPanelWidth - 8,
      kPanelHeight - kMenuHeight - 8);
    hierarchy_group_->end();
    hierarchy_group_->resizable(hierarchy_);

    attributes_group_ =
      new Fl_Group(0, kMenuHeight, kPanelWidth,
                   kPanelHeight - kMenuHeight, "Attributes");
    attributes_group_->begin();
    attributes_ = new Fl_Ged_Attributes(
      4, kMenuHeight + 4, kPanelWidth - 8,
      kPanelHeight - kMenuHeight - 8);
    attributes_group_->end();
    attributes_group_->resizable(attributes_);

    view_group_ =
      new Fl_Group(0, kMenuHeight, kPanelWidth,
                   kPanelHeight - kMenuHeight, "View");
    view_group_->begin();
    view_controls_ = new Fl_Ged_View_Controls(
      4, kMenuHeight + 4, kPanelWidth - 8,
      kPanelHeight - kMenuHeight - 8);
    view_group_->end();
    view_group_->resizable(view_controls_);

    interaction_group_ =
      new Fl_Group(0, kMenuHeight, kPanelWidth,
                   kPanelHeight - kMenuHeight, "Interact");
    interaction_group_->begin();
    interaction_ = new Fl_Ged_Interaction(
      4, kMenuHeight + 4, kPanelWidth - 8,
      kPanelHeight - kMenuHeight - 8);
    interaction_group_->end();
    interaction_group_->resizable(interaction_);
    panels_->end();

    panels_window_->end();
    panels_window_->resizable(panels_);
    panels_window_->callback([](Fl_Widget* widget, void*) { widget->hide(); });
    fltkcad::configure_accessibility(owner_);
    fltkcad::configure_accessibility(*graphics_window_);
    fltkcad::configure_accessibility(*panels_window_);

    initialize_faceplates();
    attach_views();
    terminal_->command_handler([this](const std::string& command) {
      terminal_submission_ = true;
      struct ResetSubmission {
        bool& active;
        ~ResetSubmission() { active = false; }
      } reset{terminal_submission_};
      execute_prompt(command);
    });
    terminal_->interrupt_handler([this]() {
      tcl_->interrupt();
      reset_pending_command();
      if (ged_ && ged_->dbip) execute("rtabort");
    });
    terminal_->completion_handler(
      [this](const std::string& input, std::size_t cursor) {
        return completions(input, cursor);
      });
    views_->active_handler([this](std::size_t) {
      activate_view();
      view_controls_->view(&views_->active_view());
      refresh_faceplate();
    });
    views_->changed_handler([this](std::size_t) {
      activate_view();
      view_controls_->refresh();
      refresh_faceplate();
    });
    views_->error_handler([this](std::size_t, const std::string& message) {
      report_error(message);
      terminal_->append("[display] " + message + "\n");
    });
    views_->input_handler(
      [this](std::size_t index, const Fl_Cad_View::InputEvent& event) {
        views_->active_index(index);
        activate_view();
        return interaction_->process_input(views_->view(index), event);
      });
    hierarchy_->selection_handler(
      [this](const Fl_Ged_Hierarchy::Selection& selection) {
        selected_path_ = display_path(selection.path);
        selected_object_ = selection.object;
        interaction_->selection(selected_path_, selected_object_);
        attributes_->object(selected_object_);
        refresh_faceplate();
      });
    hierarchy_->command_handler([this](const std::string& command) {
      execute(command);
    });
    hierarchy_->error_handler([this](const std::string& message) {
      report_panel_error(message);
    });
    attributes_->changed_handler([this]() {
      refresh_database_panels();
      views_->invalidate();
    });
    attributes_->error_handler([this](const std::string& message) {
      report_panel_error(message);
    });
    view_controls_->view(&views_->active_view());
    view_controls_->changed_handler([this]() { refresh_faceplate(); });
    view_controls_->error_handler([this](const std::string& message) {
      report_panel_error(message);
    });
    interaction_->command_handler(
      [this](const std::string& command) { return execute(command); });
    interaction_->selection_handler(
      [this](const Fl_Ged_Interaction::Selection& selection) {
        selected_path_ = selection.path;
        selected_object_ = selection.object;
        attributes_->object(selected_object_);
        refresh_faceplate();
      });
    interaction_->changed_handler([this]() {
      refresh_database_panels();
      views_->invalidate();
    });
    interaction_->error_handler([this](const std::string& message) {
      report_panel_error(message);
    });
    hierarchy_->ged(ged_);
    attributes_->ged(ged_);
    interaction_->ged(ged_);
    tcl_->ged(ged_);
    tcl_->changed_handler([this]() {
      refresh_database_panels();
      views_->invalidate();
      refresh_faceplate();
    });
    tcl_->output_handler([this](const std::string& text, bool) {
      flush_log();
      terminal_->append(text);
    });
    tcl_->event_processing(true);
    tcl_->ui_handler(
      [this](const std::vector<std::string>& arguments, std::string& result) {
        return handle_tcl_ui(arguments, result);
      });
    apply_preferences(load_preferences(), false);
    refresh_faceplate();
    bu_log_add_hook(capture_log, this);
    Fl::add_timeout(kProcessPollIntervalSeconds, log_timeout, this);
  }

  ~Impl()
  {
    tcl_->event_processing(false);
    Fl::remove_timeout(log_timeout, this);
    for (auto* view : views_->views()) view->faceplate_gui(nullptr);
    views_->input_handler({});
    views_->active_handler({});
    views_->changed_handler({});
    views_->error_handler({});
    terminal_->command_handler({});
    terminal_->completion_handler({});
    terminal_->interrupt_handler({});
    release_database_tools();
    close_ged();
    bu_log_delete_hook(capture_log, this);
  }

  static int capture_log(void* data, void* message)
  {
    auto& self = *static_cast<Impl*>(data);
    const char* text = static_cast<const char*>(message);
    if (!text) return 0;
    std::lock_guard<std::mutex> lock(self.log_mutex_);
    const std::size_t available = kMaxPendingLogBytes - self.pending_log_.size();
    const std::size_t length = std::char_traits<char>::length(text);
    self.pending_log_.append(text, std::min(available, length));
    self.log_truncated_ = self.log_truncated_ || length > available;
    return 0;
  }

  void flush_log()
  {
    std::string text;
    bool truncated = false;
    {
      std::lock_guard<std::mutex> lock(log_mutex_);
      text.swap(pending_log_);
      std::swap(truncated, log_truncated_);
    }
    if (!text.empty()) terminal_->append(text);
    if (truncated) terminal_->append("\n[diagnostic output truncated by terminal limit]\n");
  }

  static void log_timeout(void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.flush_log();
    const int tick = ged_edit_session_tick(self.ged_, 1);
    if (tick != 0) {
      self.views_->invalidate();
      self.refresh_faceplate();
      if (tick < 0) self.append_ged_result();
    }
    Fl::repeat_timeout(kProcessPollIntervalSeconds, log_timeout, data);
  }

  bool open_database(const std::string& path)
  {
    if (path.empty() || path.find('\0') != std::string::npos) {
      report_error("A non-empty database path is required");
      return false;
    }
    if (has_unsaved_edits()) {
      report_error("Save or discard open editor changes before replacing the database");
      terminal_->append("[error] " + last_error_ + "\n");
      return false;
    }

    struct ged* replacement = ged_open("db", path.c_str(), 1);
    if (!replacement) {
      report_error("Unable to open BRL-CAD database: " + path);
      terminal_->append("[error] " + last_error_ + "\n");
      return false;
    }

    release_database_tools();
    close_ged();
    ged_ = replacement;
    attach_views();
    bind_database();
    terminal_->append("Opened " + path + "\n");
    last_error_.clear();
    return true;
  }

  bool has_unsaved_edits() const
  {
    return ged_edit_session_mode(ged_) != GED_EDIT_VIEW ||
      (primitive_editor_ && primitive_editor_->dirty()) ||
      (combination_editor_ && combination_editor_->dirty()) ||
      (sketch_editor_ && sketch_editor_->dirty());
  }

  void release_database_tools()
  {
    primitive_editor_.reset();
    combination_editor_.reset();
    sketch_editor_.reset();
    create_primitive_.reset();
    search_.reset();
    specialist_tools_.reset();
    transform_editor_.reset();
    render_tools_.reset();
  }

  void bind_database()
  {
    reset_pending_command();
    const std::string path = ged_->dbip ? ged_->dbip->dbi_filename : "fltkged";
    owner_.copy_label(("MGED Command Window - " + path).c_str());
    graphics_window_->copy_label(("MGED Graphics Window - " + path).c_str());
    selected_object_.clear();
    selected_path_.clear();
    hierarchy_->ged(ged_);
    attributes_->ged(ged_);
    interaction_->ged(ged_);
    interaction_->selection({}, {});
    tcl_->ged(ged_);
    views_->invalidate();
    view_controls_->view(&views_->active_view());
    refresh_faceplate();
  }

  int execute(const std::string& command)
  {
    if (!ged_) {
      report_error("GED is not initialized");
      return BRLCAD_ERROR;
    }
    if (command.find('\0') != std::string::npos ||
        command.size() > terminal_->limits().max_input_bytes) {
      report_error("Command exceeds the configured input limit");
      terminal_->append("[error] " + last_error_ + "\n");
      reset_pending_command();
      return BRLCAD_ERROR;
    }

    ged_set_terminal_columns(ged_, static_cast<std::size_t>(terminal_->columns()));

    if (pending_arguments_.empty() &&
        (command == "tcl" || command.rfind("tcl ", 0) == 0)) {
      if (command == "tcl") {
        report_panel_error("Usage: tcl SCRIPT");
        return BRLCAD_ERROR;
      }
      ged_->ged_gvp = views_->active_view().view();
      const int status = tcl_->evaluate(command.substr(4));
      flush_log();
      if (tcl_->exit_requested()) return BRLCAD_OK;
      const std::string& output = tcl_->result();
      if (!output.empty()) {
        terminal_->append(output);
        if (output.back() != '\n') {
          terminal_->append("\n");
        }
      }
      if (status != 0) {
        report_error(tcl_->last_error());
      }
      return status == 0 ? BRLCAD_OK : BRLCAD_ERROR;
    }

    std::vector<std::string> current = parse_command(command);
    if (current.empty()) {
      if (pending_arguments_.empty()) {
        return BRLCAD_OK;
      }
      current.emplace_back();
    }

    std::vector<std::string> arguments = pending_arguments_;
    arguments.insert(arguments.end(), current.begin(), current.end());
    if (arguments.size() > kMaxCommandArguments ||
        argument_bytes(arguments) > kMaxPendingCommandBytes) {
      report_error("Incremental command exceeds the application limit");
      terminal_->append("[error] " + last_error_ + "\n");
      reset_pending_command();
      return BRLCAD_ERROR;
    }

    std::vector<const char*> argv;
    argv.reserve(arguments.size() + 1);
    for (const auto& argument : arguments) {
      argv.push_back(argument.c_str());
    }
    argv.push_back(nullptr);

    bu_vls_trunc(ged_->ged_result_str, 0);
    ged_->ged_gvp = views_->active_view().view();
    const int result = ged_edit_session_exec(ged_, static_cast<int>(arguments.size()), argv.data());
    flush_log();
    const bool needs_more = (result & GED_MORE) != 0;
    if (needs_more) {
      begin_pending_history();
      pending_arguments_ = std::move(arguments);
      set_continuation_prompt();
    } else {
      append_ged_result();
      finish_pending_history();
      pending_arguments_.clear();
      terminal_->prompt(kDefaultPrompt);
    }
    bu_vls_trunc(ged_->ged_result_str, 0);
    views_->invalidate();
    refresh_faceplate();
    if (!needs_more && result == BRLCAD_OK && may_change_database(arguments)) {
      refresh_database_panels();
    }
    return result & ~GED_EXIT;
  }

  int execute_prompt(const std::string& command)
  {
    if (command_language_ == FltkGedWindow::CommandLanguage::Ged)
      return execute(command);
    if (pending_script_.empty() &&
        (command == "tcl" || command.rfind("tcl ", 0) == 0))
      return execute(command);
    if (!tcl_->available()) {
      report_panel_error(tcl_->last_error());
      return BRLCAD_ERROR;
    }
    if (command.find('\0') != std::string::npos ||
        pending_script_.size() + command.size() > terminal_->limits().max_input_bytes) {
      pending_script_.clear();
      reset_pending_command();
      report_panel_error("Tcl command exceeds the configured input limit");
      return BRLCAD_ERROR;
    }
    if (!pending_script_.empty()) pending_script_ += pending_tcl_more_ ? ' ' : '\n';
    pending_script_ += command;
    if (!tcl_->command_complete(pending_script_)) {
      begin_pending_history();
      terminal_->prompt("? ");
      return GED_MORE;
    }
    ged_set_terminal_columns(ged_, static_cast<std::size_t>(terminal_->columns()));
    ged_->ged_gvp = views_->active_view().view();
    const std::string script = std::exchange(pending_script_, {});
    const int status = tcl_->evaluate(script);
    flush_log();
    if (status != 0 && (tcl_->command_status() & GED_MORE)) {
      pending_script_ = script;
      pending_tcl_more_ = true;
      begin_pending_history();
      terminal_->prompt(tcl_->result().empty() ? kContinuationPrompt : tcl_->result());
      return GED_MORE;
    }
    pending_tcl_more_ = false;
    finish_pending_history(&script);
    if (!tcl_->exit_requested()) {
      const std::string& output = status == 0 ? tcl_->result() : tcl_->last_error();
      if (!output.empty()) terminal_->append(output + (output.back() == '\n' ? "" : "\n"));
    }
    terminal_->prompt(kDefaultPrompt);
    views_->invalidate();
    refresh_faceplate();
    refresh_database_panels();
    return status == 0 || tcl_->exit_requested() ? BRLCAD_OK : BRLCAD_ERROR;
  }

  void command_language(FltkGedWindow::CommandLanguage language)
  {
    if (language == FltkGedWindow::CommandLanguage::Tcl && !tcl_->available()) {
      report_panel_error(tcl_->last_error());
      return;
    }
    pending_script_.clear();
    reset_pending_command();
    command_language_ = language;
    if (const Fl_Menu_Item* item = menu_->find_item(
          FLTKGED_TR("&Settings/Command Language/&Tcl"))) {
      if (language == FltkGedWindow::CommandLanguage::Tcl)
        const_cast<Fl_Menu_Item*>(item)->setonly(menu_->menu());
    }
    if (const Fl_Menu_Item* item = menu_->find_item(
          FLTKGED_TR("&Settings/Command Language/&Raw libged"))) {
      if (language == FltkGedWindow::CommandLanguage::Ged)
        const_cast<Fl_Menu_Item*>(item)->setonly(menu_->menu());
    }
    sync_menus();
  }

  static void tcl_prompt_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->command_language(FltkGedWindow::CommandLanguage::Tcl);
  }

  static void ged_prompt_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->command_language(FltkGedWindow::CommandLanguage::Ged);
  }

  Fl_Cad_Terminal::Completion completions(const std::string& input,
                                          std::size_t cursor)
  {
    Fl_Cad_Terminal::Completion result;
    cursor = std::min(cursor, input.size());
    std::size_t token_begin = cursor;
    while (token_begin > 0 && !is_space(input[token_begin - 1])) {
      --token_begin;
    }
    const std::string seed = input.substr(token_begin, cursor - token_begin);
    const bool command_token = std::all_of(
      input.begin(), input.begin() + static_cast<std::ptrdiff_t>(token_begin),
      is_space);

    const char** matches = nullptr;
    int count = 0;
    result.begin = token_begin;
    result.end = cursor;
    if (command_token) {
      count = ged_cmd_completions(&matches, seed.c_str());
    } else if (ged_ && ged_->dbip) {
      struct bu_vls prefix = BU_VLS_INIT_ZERO;
      count = ged_geom_completions(&matches, &prefix, ged_->dbip, seed.c_str());
      const std::size_t prefix_size = bu_vls_strlen(&prefix);
      if (prefix_size <= cursor) {
        result.begin = cursor - prefix_size;
      }
      bu_vls_free(&prefix);
    }

    if (count > 0 && matches) {
      const std::size_t limit = terminal_->limits().max_completion_candidates;
      const int accepted = std::min(
        count, static_cast<int>(std::min(
                 limit, static_cast<std::size_t>(std::numeric_limits<int>::max()))));
      result.candidates.reserve(static_cast<std::size_t>(accepted));
      for (int index = 0; index < accepted; ++index) {
        result.candidates.emplace_back(matches[index]);
      }
      bu_argv_free(count, const_cast<char**>(matches));
      if (accepted != count) {
        report_error("Completion results exceeded the configured limit");
      }
    }
    return result;
  }

  void initialize_faceplates()
  {
    for (std::size_t i = 0; i < faceplates_.size(); ++i) {
      auto& faceplate = faceplates_[i];
      faceplate.owner = this;
      faceplate.widget = &views_->view(i);
      dm_faceplate_controls_init(&faceplate.menu, DM_FACEPLATE_MENU);
      faceplate.menu.data = &faceplate;
      faceplate.menu.item = faceplate_item;
      faceplate.menu.select = faceplate_select;
      dm_faceplate_controls_init(&faceplate.sliders, DM_FACEPLATE_SLIDER);
      faceplate.sliders.data = &faceplate;
      faceplate.sliders.item = slider_item;
      faceplate.sliders.select = slider_select;
      faceplate.gui.menu = &faceplate.menu;
      faceplate.widget->view()->gv_s->gv_view_params.draw = 1;
    }
  }

  static int faceplate_item(void* data, std::size_t row, struct dm_faceplate_item* item)
  {
    auto& self = *static_cast<Faceplate*>(data)->owner;
    if (ged_edit_session_mode(self.ged_) != GED_EDIT_VIEW) {
      if (row && row <= std::size(kEditMenuActions)) {
        item->label = kEditMenuActions[row - 1].label;
      } else {
        const std::size_t menu_row = row ? row - std::size(kEditMenuActions) : 0;
        item->label = ged_edit_session_menu_label(self.ged_, menu_row);
        if (!item->label && !row) item->label = "Solid Edit";
      }
      if (!item->label) return 0;
      VSET(item->color, 255, 255, row ? 0 : 255);
      return 1;
    }
    if (row >= std::size(kFaceplateMenu)) return 0;
    item->label = kFaceplateMenu[row].label;
    VSET(item->color, 255, 255, row ? 0 : 255);
    return 1;
  }

  static void faceplate_select(void* data, std::size_t row, fastf_t, int activate)
  {
    if (!activate || !row) return;
    auto& faceplate = *static_cast<Faceplate*>(data);
    auto& self = *faceplate.owner;
    self.activate_faceplate(faceplate);
    if (ged_edit_session_mode(self.ged_) != GED_EDIT_VIEW) {
      if (row <= std::size(kEditMenuActions)) {
        self.execute(kEditMenuActions[row - 1].command);
      } else if (const char* label = ged_edit_session_menu_label(
                   self.ged_, row - std::size(kEditMenuActions))) {
        self.execute("press " + quote(label));
      }
      return;
    }
    if (row >= std::size(kFaceplateMenu)) return;
    auto* view = faceplate.widget->view();
    vect_t aet;
    switch (kFaceplateMenu[row].action) {
      case FaceplateAction::Title: return;
      case FaceplateAction::Oblique: VSET(aet, 35, 25, 0); break;
      case FaceplateAction::Top: VSET(aet, 270, 90, 0); break;
      case FaceplateAction::Right: VSET(aet, 270, 0, 0); break;
      case FaceplateAction::Front: VSET(aet, 0, 0, 0); break;
      case FaceplateAction::Isometric: VSET(aet, 45, 45, 0); break;
      case FaceplateAction::Restore: self.step_stored_view(0); return;
      case FaceplateAction::Save: self.add_stored_view(); return;
      case FaceplateAction::Adc:
        view->gv_s->gv_adc.draw = !view->gv_s->gv_adc.draw;
        return;
      case FaceplateAction::Sliders:
        faceplate.sliders_enabled = !faceplate.sliders_enabled;
        self.refresh_faceplate();
        return;
      case FaceplateAction::ZoomIn: bv_view_set_scale(view, view->gv_scale * 0.5); return;
      case FaceplateAction::ZoomOut: bv_view_set_scale(view, view->gv_scale * 2.0); return;
      default: return;
    }
    bv_view_set_aet(view, aet);
    faceplate.values.fill(0);
  }

  static int slider_item(void* data, std::size_t row, struct dm_faceplate_item* item)
  {
    auto& faceplate = *static_cast<Faceplate*>(data);
    if (row >= std::size(kSliderLabels)) return 0;
    item->label = kSliderLabels[row];
    item->value = faceplate.values[row];
    VSET(item->color, 255, 255, 0);
    return 1;
  }

  void activate_faceplate(const Faceplate& faceplate)
  {
    for (std::size_t i = 0; i < faceplates_.size(); ++i) {
      if (&faceplates_[i] == &faceplate) {
        views_->active_index(i);
        activate_view();
        return;
      }
    }
  }

  static void slider_select(void* data, std::size_t row, fastf_t value, int activate)
  {
    auto& faceplate = *static_cast<Faceplate*>(data);
    if (!activate || row >= faceplate.values.size()) return;
    faceplate.owner->activate_faceplate(faceplate);
    const double dead_zone = faceplate.sliders.dead_zone;
    value = value > dead_zone ? value - dead_zone :
      value < -dead_zone ? value + dead_zone : 0;
    const double delta = value - faceplate.values[row];
    faceplate.values[row] = value;
    auto* view = faceplate.widget->view();
    constexpr std::size_t scale_row = Z + 1;
    constexpr double scale_range = 4.0;
    constexpr double rotation_range = 180.0;
    if (row < scale_row) {
      vect_t offset = VINIT_ZERO;
      offset[row] = -delta;
      vect_t model_offset;
      point_t center;
      MAT4X3VEC(model_offset, view->gv_view2model, offset);
      MAT_DELTAS_GET_NEG(center, view->gv_center);
      VADD2(center, center, model_offset);
      bv_view_set_center_vec(view, center);
    } else if (row == scale_row) {
      bv_view_set_scale(view, view->gv_scale * std::pow(scale_range, -delta));
    } else {
      vect_t angles = VINIT_ZERO;
      angles[row - scale_row - 1] = delta * rotation_range;
      mat_t rotation, result;
      bn_mat_angles(rotation, V3ARGS(angles));
      bn_mat_mul(result, rotation, view->gv_rotation);
      bv_view_set_rotation(view, result);
    }
  }

  static void faceplate_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    auto& view = self.views_->active_view();
    auto& params = view.view()->gv_s->gv_view_params;
    params.draw = !params.draw;
    view.invalidate();
    self.refresh_faceplate();
  }

  static void faceplate_gui_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    auto& faceplate = self.faceplates_[self.views_->active_index()];
    faceplate.enabled = !faceplate.enabled;
    if (faceplate.enabled) faceplate.widget->view()->gv_s->gv_view_params.draw = 1;
    self.refresh_faceplate();
  }

  void attach_views()
  {
    if (!ged_ || !views_) {
      return;
    }
    for (auto* widget : views_->views()) {
      if (!widget || !widget->view()) {
        continue;
      }
      widget->ged(ged_);
      bv_set_add_view(&ged_->ged_views, widget->view());
      if (ged_->dbip) {
        widget->view()->gv_base2local = ged_->dbip->dbi_base2local;
        widget->view()->gv_local2base = ged_->dbip->dbi_local2base;
      }
    }
    ged_->ged_gvp = views_->active_view().view();
    ged_->vlist_ctx = this;
    ged_->ged_create_vlist_scene_obj_callback = activate_scene_object;
    ged_->ged_create_io_handler = create_io_handler;
    ged_->ged_delete_io_handler = delete_io_handler;
    ged_->ged_io_data = this;
    ged_clbk_set(ged_, "clear", BU_CLBK_POST, clear_terminal, this);
    // libged keys callbacks by implementation, so this also handles q and exit.
    ged_clbk_set(ged_, "quit", BU_CLBK_POST, quit_command, this);
    ged_clbk_set(ged_, "opendb", BU_CLBK_PRE, before_open_database, this);
    ged_clbk_set(ged_, "opendb", BU_CLBK_POST, after_open_database, this);
    ged_clbk_set(ged_, "closedb", BU_CLBK_PRE, before_close_database, this);
    ged_clbk_set(ged_, "closedb", BU_CLBK_POST, after_close_database, this);
    ged_clbk_set(ged_, "search", BU_CLBK_DURING, search_execute, this);
    for (const char* command : {"draw", "E", "ev"}) {
      ged_clbk_set(ged_, command, BU_CLBK_PRE, before_draw, this);
      ged_clbk_set(ged_, command, BU_CLBK_POST, after_draw, this);
    }
    if (ged_->ged_fbs) {
      ged_->ged_fbs->fbs_open_ipc_client_handler = open_fb_client;
      ged_->ged_fbs->fbs_close_ipc_client_handler = close_fb_client;
      ged_->ged_fbs->fbs_close_server_handler = close_fb_server;
      ged_->ged_fbs->fbs_callback = framebuffer_changed;
      ged_->ged_fbs->fbs_clientData = this;
    }
  }

  void activate_view()
  {
    if (ged_ && views_) {
      ged_->ged_gvp = views_->active_view().view();
    }
  }

  static int clear_terminal(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (ged_results_ret(self.ged_->ged_results) != BRLCAD_OK) return BRLCAD_OK;
    self.flush_log();
    self.terminal_->clear();
    return BRLCAD_OK;
  }

  static int quit_command(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (ged_results_ret(self.ged_->ged_results) & GED_EXIT) self.request_quit();
    return BRLCAD_OK;
  }

  static int before_open_database(int argc, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.database_before_open_ = self.ged_->dbip;
    if (argc <= 1) return BRLCAD_OK;
    return before_close_database(0, nullptr, nullptr, data);
  }

  static int before_close_database(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.has_unsaved_edits()) return BRLCAD_OK;
    bu_vls_strcat(self.ged_->ged_result_str,
      "Save or discard open editor changes before replacing or closing the database\n");
    return BRLCAD_ERROR;
  }

  static int after_open_database(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.ged_->dbip && self.ged_->dbip != self.database_before_open_) {
      for (auto* widget : self.views_->views()) {
        widget->view()->gv_base2local = self.ged_->dbip->dbi_base2local;
        widget->view()->gv_local2base = self.ged_->dbip->dbi_local2base;
      }
      self.bind_database();
    }
    return BRLCAD_OK;
  }

  static int after_close_database(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_->dbip) {
      self.release_database_tools();
      self.bind_database();
    }
    return BRLCAD_OK;
  }

  static int search_execute(int argc, const char** argv, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const std::string output = bu_vls_cstr(self.ged_->ged_result_str);
    bu_vls_trunc(self.ged_->ged_result_str, 0);
    char** arguments = bu_argv_dup(argc, argv);
    const int status = ged_edit_session_exec(self.ged_, argc, const_cast<const char**>(arguments));
    bu_argv_free(argc, arguments);
    self.flush_log();
    self.append_ged_result();
    bu_vls_sprintf(self.ged_->ged_result_str, "%s", output.c_str());
    return (status & ~GED_EXIT) == BRLCAD_OK;
  }

  static int before_draw(int argc, const char** argv, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    bool keep_view = false;
    for (int index = 1; index < argc; ++index) {
      const char* option = argv[index];
      if (option[0] != '-' || BU_STR_EQUAL(option, "--")) break;
      for (++option; *option; ++option) {
        if (*option == 'R') keep_view = true;
        // These legacy draw options consume the rest of this token or the
        // next argument.  An R inside their value is not a resize flag.
        constexpr const char* options_with_values = "mxCPL";
        if (std::strchr(options_with_values, *option)) {
          if (!option[1]) ++index;
          break;
        }
      }
    }
    self.draw_autoview_.push_back(!keep_view && ged_who_argc(self.ged_) == 0);
    return BRLCAD_OK;
  }

  static int after_draw(int, const char**, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.draw_autoview_.empty()) return BRLCAD_OK;
    const bool resize = self.draw_autoview_.back();
    self.draw_autoview_.pop_back();
    if (!resize || ged_results_ret(self.ged_->ged_results) != BRLCAD_OK ||
        ged_who_argc(self.ged_) == 0) {
      return BRLCAD_OK;
    }

    // MGED fits a newly populated display, preserving later camera changes.
    // Keep command diagnostics intact while framing the shared geometry.
    const std::string output = bu_vls_cstr(self.ged_->ged_result_str);
    const char* arguments[] = {"autoview"};
    const int status = ged_exec_autoview(self.ged_, 1, arguments);
    if (status != BRLCAD_OK) {
      const std::string error = bu_vls_cstr(self.ged_->ged_result_str);
      bu_vls_sprintf(self.ged_->ged_result_str, "%s%s", output.c_str(), error.c_str());
      return status;
    }
    bu_vls_sprintf(self.ged_->ged_result_str, "%s", output.c_str());
    const auto* active = self.ged_->ged_gvp;
    point_t center;
    bv_view_get_center_vec(active, center);
    for (auto* widget : self.views_->views()) {
      if (widget->view() == active) continue;
      bv_view_set_center_vec(widget->view(), center);
      bv_view_set_scale(widget->view(), bv_view_get_scale(active));
      bv_update(widget->view());
    }
    return BRLCAD_OK;
  }

  int handle_tcl_ui(const std::vector<std::string>& arguments,
                    std::string& result)
  {
    const auto usage = [&result]() {
      result = "usage: fltkged database | selection ?object? | view ACTION | show PANEL";
      return 1;
    };
    if (arguments.empty()) {
      return usage();
    }
    if (arguments[0] == "database" && arguments.size() == 1) {
      result = ged_ && ged_->dbip && ged_->dbip->dbi_filename
                 ? ged_->dbip->dbi_filename
                 : "";
      return 0;
    }
    if (arguments[0] == "selection") {
      if (arguments.size() == 1) {
        result = selected_object_;
        return 0;
      }
      if (arguments.size() != 2 || !ged_ || !ged_->dbip ||
          !db_lookup(ged_->dbip, arguments[1].c_str(), LOOKUP_QUIET)) {
        result = "selection must name an object in the open database";
        return 1;
      }
      selected_object_ = arguments[1];
      selected_path_ = "/" + selected_object_;
      attributes_->object(selected_object_);
      interaction_->selection(selected_path_, selected_object_);
      refresh_faceplate();
      result = selected_object_;
      return 0;
    }
    if (arguments[0] == "view" && arguments.size() == 2) {
      if (arguments[1] == "redraw") {
        views_->invalidate();
      } else if (arguments[1] == "autoview") {
        execute("autoview");
      } else if (arguments[1] == "single") {
        views_->layout(Fl_Cad_Quad_View::Layout::Single);
      } else if (arguments[1] == "quad") {
        views_->layout(Fl_Cad_Quad_View::Layout::Quad);
      } else {
        result = "view action must be redraw, autoview, single, or quad";
        return 1;
      }
      result = arguments[1];
      return 0;
    }
    if (arguments[0] == "show" && arguments.size() == 2) {
      const std::string& panel = arguments[1];
      if (panel == "preferences") {
        preferences_menu(nullptr, this);
      } else if (panel == "primitive") {
        edit_primitive_menu(nullptr, this);
      } else if (panel == "combination") {
        edit_combination_menu(nullptr, this);
      } else if (panel == "sketch") {
        edit_sketch_menu(nullptr, this);
      } else if (panel == "render") {
        render_tools_menu(nullptr, this);
      } else if (panel == "specialist") {
        specialist_tools_menu(nullptr, this);
      } else if (panel == "help") {
        documentation_menu(nullptr, this);
      } else if (panel == "about") {
        about_menu(nullptr, this);
      } else {
        result = "unknown FLTK panel: " + panel;
        return 1;
      }
      result = panel;
      return 0;
    }
    return usage();
  }

  static void activate_scene_object(void*, struct bv_scene_obj* object)
  {
    if (object) {
      object->s_flag = UP;
    }
  }

  void close_ged()
  {
    if (!ged_) {
      return;
    }
    if (tcl_) {
      tcl_->ged(nullptr);
    }
    ged_subprocesses_terminate(ged_);
    clear_io_listeners();
    close_framebuffer_server();
    // Release the legacy display lists before clearing view-local objects.
    const char* arguments[] = {"zap"};
    if (ged_exec_zap(ged_, 1, arguments) != BRLCAD_OK) {
      bu_log("Unable to clear geometry during shutdown: %s\n",
             bu_vls_cstr(ged_->ged_result_str));
    }
    if (views_) {
      for (auto* widget : views_->views()) {
        if (!widget || !widget->view()) {
          continue;
        }
        bv_clear(widget->view(), 0);
        widget->ged(nullptr);
        bv_set_rm_view(&ged_->ged_views, widget->view());
        if (ged_->ged_gvp == widget->view()) {
          ged_->ged_gvp = nullptr;
        }
      }
    }
    ged_close(ged_);
    ged_ = nullptr;
  }

  void append_ged_result()
  {
    const std::size_t size = bu_vls_strlen(ged_->ged_result_str);
    if (size == 0) {
      return;
    }
    const std::size_t limit = terminal_->limits().max_append_bytes;
    const std::size_t accepted = std::min(size, limit);
    std::string output(bu_vls_cstr(ged_->ged_result_str), accepted);
    terminal_->append(output);
    if (!output.empty() && output.back() != '\n') {
      terminal_->append("\n");
    }
    if (accepted != size) {
      terminal_->append(kOutputTruncation);
      report_error("GED output exceeded the configured terminal limit");
    }
  }

  void set_continuation_prompt()
  {
    const std::size_t size = bu_vls_strlen(ged_->ged_result_str);
    if (size == 0) {
      terminal_->prompt(kContinuationPrompt);
      return;
    }
    const std::size_t accepted = std::min(size, kMaxPromptBytes);
    terminal_->prompt(std::string(bu_vls_cstr(ged_->ged_result_str), accepted));
    if (accepted != size) {
      report_error("GED continuation prompt exceeded the application limit");
    }
  }

  void finish_pending_history(const std::string* script = nullptr)
  {
    if (!pending_history_begin_) {
      return;
    }
    const std::size_t history_size = terminal_->history_size();
    if (history_size > 0) {
      if (script) terminal_->consolidate_history(*pending_history_begin_, history_size, *script);
      else terminal_->consolidate_history(*pending_history_begin_, history_size);
    }
    pending_history_begin_.reset();
  }

  void begin_pending_history()
  {
    if (terminal_submission_ && !pending_history_begin_ && terminal_->history_size())
      pending_history_begin_ = terminal_->history_size() - 1;
  }

  void reset_pending_command()
  {
    pending_script_.clear();
    pending_tcl_more_ = false;
    pending_arguments_.clear();
    pending_history_begin_.reset();
    terminal_->prompt(kDefaultPrompt);
  }

  void report_error(std::string message)
  {
    last_error_ = std::move(message);
  }

  void report_panel_error(const std::string& message)
  {
    report_error(message);
    terminal_->append("[error] " + message + "\n");
  }

  void refresh_database_panels()
  {
    hierarchy_->refresh();
    if (!selected_object_.empty() && ged_ && ged_->dbip &&
        db_lookup(ged_->dbip, selected_object_.c_str(), LOOKUP_QUIET)) {
      attributes_->object(selected_object_);
      attributes_->refresh();
    } else {
      selected_object_.clear();
      selected_path_.clear();
      interaction_->selection({}, {});
      attributes_->object({});
    }
  }

  void refresh_faceplate()
  {
    for (auto& faceplate : faceplates_) {
      faceplate.gui.sliders = faceplate.sliders_enabled ? &faceplate.sliders : nullptr;
      const bool visible = faceplate.enabled &&
        faceplate.widget->view()->gv_s->gv_view_params.draw;
      faceplate.widget->faceplate_gui(visible ? &faceplate.gui : nullptr);
    }
    sync_menus();
  }

  void add_stored_view()
  {
    const std::size_t view_index = views_->active_index();
    auto& ring = stored_views_[view_index];
    if (ring.size() == kMaxStoredViews) {
      report_panel_error("The active view already has the maximum of ten stored views");
      return;
    }
    struct bview* view = views_->active_view().view();
    if (!view) {
      report_panel_error("No active view is available");
      return;
    }
    ViewSnapshot snapshot;
    bv_view_get_center_vec(view, snapshot.center);
    bv_view_get_aet(view, snapshot.aet);
    snapshot.size = bv_view_get_size(view);
    snapshot.perspective = bv_view_get_perspective(view);
    ring.push_back(snapshot);
    stored_view_indices_[view_index] = ring.size() - 1;
    refresh_faceplate();
  }

  void step_stored_view(int direction)
  {
    const std::size_t view_index = views_->active_index();
    auto& ring = stored_views_[view_index];
    if (ring.empty()) {
      report_panel_error("The active view has no stored views");
      return;
    }
    std::size_t& current = stored_view_indices_[view_index];
    if (current >= ring.size()) {
      current = 0;
    } else if (direction > 0) {
      current = (current + 1) % ring.size();
    } else if (direction < 0) {
      current = (current + ring.size() - 1) % ring.size();
    }
    restore_stored_view(ring[current]);
  }

  void delete_stored_view()
  {
    const std::size_t view_index = views_->active_index();
    auto& ring = stored_views_[view_index];
    std::size_t& current = stored_view_indices_[view_index];
    if (ring.empty() || current >= ring.size()) {
      report_panel_error("The active view has no current stored view");
      return;
    }
    ring.erase(ring.begin() + static_cast<std::ptrdiff_t>(current));
    if (ring.empty()) {
      current = std::numeric_limits<std::size_t>::max();
    } else {
      current %= ring.size();
      restore_stored_view(ring[current]);
    }
    refresh_faceplate();
  }

  void restore_stored_view(const ViewSnapshot& snapshot)
  {
    struct bview* view = views_->active_view().view();
    if (!view) {
      return;
    }
    bv_view_set_center_vec(view, snapshot.center);
    bv_view_set_aet(view, snapshot.aet);
    bv_view_set_size(view, snapshot.size);
    bv_view_set_perspective(view, snapshot.perspective);
    bv_update(view);
    view_controls_->refresh();
    views_->active_view().invalidate();
    refresh_faceplate();
  }

  Fl_Ged_Preferences::Settings load_preferences() const
  {
    Fl_Ged_Preferences::Settings settings;
    Fl_Preferences preferences(Fl_Preferences::USER_L, kPreferencesVendor,
                               kPreferencesApplication);
    const auto load_color = [&preferences](const char* key,
                                           Fl_Cad_View::Color fallback) {
      int r = fallback.red, g = fallback.green, b = fallback.blue;
      preferences.get((std::string(key) + "_red").c_str(), r, r);
      preferences.get((std::string(key) + "_green").c_str(), g, g);
      preferences.get((std::string(key) + "_blue").c_str(), b, b);
      return preference_color(r, g, b);
    };
    settings.background_first = load_color("background_first", settings.background_first);
    settings.background_second = load_color("background_second", settings.background_second);
    settings.terminal_foreground = load_color("terminal_foreground", settings.terminal_foreground);
    settings.terminal_background = load_color("terminal_background", settings.terminal_background);
    preferences.get("faceplate_font_size", settings.faceplate_font_size, 0);
    settings.faceplate_font_size = std::clamp(
      settings.faceplate_font_size, 0, Fl_Ged_Preferences::maximum_font_size);
    int history = static_cast<int>(settings.command_history);
    int scrollback = settings.scrollback_lines;
    preferences.get("command_history", history, history);
    preferences.get("scrollback_lines", scrollback, scrollback);
    int wrap_history = settings.wrap_history;
    preferences.get("wrap_history", wrap_history, wrap_history);
    settings.wrap_history = wrap_history != 0;
    settings.command_history = static_cast<std::size_t>(
      std::clamp(history, 1, kMaximumPreferenceHistory));
    settings.scrollback_lines =
      std::clamp(scrollback, 1, kMaximumPreferenceScrollback);
    return settings;
  }

  static Fl_Cad_View::Color preference_color(int red, int green, int blue)
  {
    constexpr int maximum = std::numeric_limits<unsigned char>::max();
    return {static_cast<unsigned char>(std::clamp(red, 0, maximum)),
            static_cast<unsigned char>(std::clamp(green, 0, maximum)),
            static_cast<unsigned char>(std::clamp(blue, 0, maximum))};
  }

  void apply_preferences(const Fl_Ged_Preferences::Settings& settings,
                         bool save)
  {
    preferences_settings_ = settings;
    for (auto* view : views_->views()) {
      view->background(settings.background_first, settings.background_second);
      view->view()->gv_s->gv_view_params.font_size = settings.faceplate_font_size;
    }
    const auto fl_color = [](const Fl_Cad_View::Color& color) {
      return fl_rgb_color(color.red, color.green, color.blue);
    };
    terminal_->colors(fl_color(settings.terminal_foreground),
                       fl_color(settings.terminal_background));
    auto terminal_limits = terminal_->limits();
    terminal_limits.command_history = settings.command_history;
    terminal_limits.scrollback_lines = settings.scrollback_lines;
    terminal_->limits(terminal_limits);
    terminal_->wrap_history(settings.wrap_history);

    if (!save) {
      return;
    }
    Fl_Preferences preferences(Fl_Preferences::USER_L, kPreferencesVendor,
                               kPreferencesApplication);
    bool saved = true;
    const auto save_value = [&preferences, &saved](const char* key, int value) {
      if (!preferences.set(key, value)) {
        saved = false;
      }
    };
    const auto save_color = [&save_value](const char* key, const Fl_Cad_View::Color& color) {
      save_value((std::string(key) + "_red").c_str(), color.red);
      save_value((std::string(key) + "_green").c_str(), color.green);
      save_value((std::string(key) + "_blue").c_str(), color.blue);
    };
    save_color("background_first", settings.background_first);
    save_color("background_second", settings.background_second);
    save_color("terminal_foreground", settings.terminal_foreground);
    save_color("terminal_background", settings.terminal_background);
    save_value("faceplate_font_size", settings.faceplate_font_size);
    save_value("command_history", static_cast<int>(settings.command_history));
    save_value("scrollback_lines", settings.scrollback_lines);
    save_value("wrap_history", settings.wrap_history);
    if (saved && preferences.flush() < 0) {
      saved = false;
    }
    if (!saved) {
      report_panel_error("Unable to save one or more application preferences");
    }
  }

  void resize(int width, int height)
  {
    menu_->resize(0, 0, width, kMenuHeight);
    const int content_height = std::max(1, height - kMenuHeight);
    main_tile_->resize(0, kMenuHeight, width, content_height);
  }

  void sync_menus()
  {
    const int index = menu_->find_index(FLTKGED_TR("&Settings/&Combined Windows"));
    if (index >= 0) {
      if (combined_) menu_->mode(index, FL_MENU_TOGGLE | FL_MENU_VALUE);
      else menu_->mode(index, FL_MENU_TOGGLE);
    }
    if (faceplates_[0].widget) {
      const auto set_toggle = [this](const char* name, bool enabled) {
        const int item = menu_->find_index(name);
        if (item >= 0) menu_->mode(item, FL_MENU_TOGGLE | (enabled ? FL_MENU_VALUE : 0));
      };
      set_toggle(FLTKGED_TR("&Modes/&Faceplate"),
                 views_->active_view().view()->gv_s->gv_view_params.draw);
      set_toggle(FLTKGED_TR("&Modes/Faceplate &GUI"),
                 faceplates_[views_->active_index()].enabled);
    }
    graphics_menu_->copy(menu_->menu());
  }

  void combined(bool value)
  {
    if (combined_ == value) return;
    views_->parent()->remove(views_);
    combined_ = value;
    if (value) {
      graphics_window_->hide();
      main_tile_->insert(*views_, 0);
      const int content_height = main_tile_->h();
      const int view_height = content_height * kInitialViewPercent / 100;
      views_->resize(0, kMenuHeight, owner_.w(), view_height);
      terminal_->resize(0, kMenuHeight + view_height, owner_.w(),
                          content_height - view_height);
      main_tile_->resizable(views_);
    } else {
      graphics_window_->add(views_);
      views_->resize(0, kMenuHeight, graphics_window_->w(),
                       graphics_window_->h() - kMenuHeight);
      graphics_window_->resizable(views_);
      terminal_->resize(0, kMenuHeight, owner_.w(), main_tile_->h());
      main_tile_->resizable(terminal_);
      if (owner_.shown()) graphics_window_->show();
    }
    main_tile_->init_sizes();
    graphics_window_->init_sizes();
    sync_menus();
    owner_.redraw();
  }

  void hide_windows()
  {
    graphics_window_->hide();
    panels_window_->hide();
    if (help_) help_->hide();
    if (about_) about_->hide();
    if (primitive_editor_) primitive_editor_->hide();
    if (combination_editor_) combination_editor_->hide();
    if (sketch_editor_) sketch_editor_->hide();
    if (create_primitive_) create_primitive_->hide();
    if (search_) search_->hide();
    if (render_tools_) render_tools_->hide();
    if (transform_editor_) transform_editor_->hide();
    if (specialist_tools_) specialist_tools_->hide();
    if (preferences_editor_) preferences_editor_->hide();
  }

  static void combined_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.combined(!self.combined_);
  }

  static void panels_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->panels_window_->show();
  }

  static void command_window_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.owner_.show();
    self.terminal_->output().take_focus();
  }

  static void graphics_window_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.combined_) self.graphics_window_->show();
  }

  static void navigate_mode_menu(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->interaction_->mode(Fl_Ged_Interaction::Mode::Navigate); }
  static void select_mode_menu(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->interaction_->mode(Fl_Ged_Interaction::Mode::Select); }
  static void query_mode_menu(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->interaction_->mode(Fl_Ged_Interaction::Mode::QueryRay); }
  static void zoom_mode_menu(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->interaction_->mode(Fl_Ged_Interaction::Mode::RectangleZoom); }
  static void rectangle_mode_menu(Fl_Widget*, void* data)
  { static_cast<Impl*>(data)->interaction_->mode(Fl_Ged_Interaction::Mode::RectangleSelect); }

  bool request_quit()
  {
    if (!has_unsaved_edits()) {
      owner_.hide();
      return true;
    }
    const int decision = fl_choice(
      "Save pending geometry edits before quitting?", "Cancel", "Discard",
      "Save");
    if (decision == 0) {
      return false;
    }
    if (decision == 2) {
      if ((primitive_editor_ && primitive_editor_->dirty() &&
           !primitive_editor_->save()) ||
          (combination_editor_ && combination_editor_->dirty() &&
           !combination_editor_->save()) ||
          (sketch_editor_ && sketch_editor_->dirty() &&
           !sketch_editor_->save()) ||
          (ged_edit_session_mode(ged_) != GED_EDIT_VIEW &&
           execute("accept") != BRLCAD_OK)) {
        report_panel_error("An editor could not save its changes; quit was cancelled");
        fl_alert("An editor could not save its changes. Quit was cancelled.");
        return false;
      }
    }
    ged_edit_session_end(ged_);
    tcl_->interrupt();
    owner_.hide();
    return true;
  }

  static void create_io_handler(struct ged_subprocess* process,
                                bu_process_io_t type,
                                ged_io_func_t callback, void* callback_data)
  {
    if (!process || !process->p || !process->gedp ||
        !process->gedp->ged_io_data || !callback) {
      return;
    }
    auto& self = *static_cast<Impl*>(process->gedp->ged_io_data);
    self.add_io_listener(process, type, callback, callback_data);
  }

  static void delete_io_handler(struct ged_subprocess* process,
                                bu_process_io_t type)
  {
    if (!process || !process->gedp || !process->gedp->ged_io_data) {
      return;
    }
    auto& self = *static_cast<Impl*>(process->gedp->ged_io_data);
    self.remove_io_listener(process, type);
  }

  void add_io_listener(struct ged_subprocess* process, bu_process_io_t type,
                       ged_io_func_t callback, void* callback_data)
  {
    const int descriptor = bu_process_fileno(process->p, type);
    if (descriptor < 0) {
      report_error("Unable to monitor a GED subprocess stream");
      return;
    }

    const IoKey key{process, static_cast<int>(type)};
    auto old = io_listeners_.find(key);
    if (old != io_listeners_.end()) {
      stop_io_listener(*old->second);
      io_listeners_.erase(old);
    }

    auto listener = std::make_unique<IoListener>();
    listener->owner = this;
    listener->process = process;
    listener->type = type;
    listener->callback = callback;
    listener->callback_data = callback_data;
    listener->descriptor = descriptor;
    listener->polling = process_polling_enabled();
    if (type != BU_PROCESS_STDIN) {
      listener->watching = true;
      if (listener->polling) {
        Fl::add_timeout(kProcessPollIntervalSeconds, io_poll, listener.get());
      } else {
        Fl::add_fd(descriptor, FL_READ, io_ready, listener.get());
      }
    }
    io_listeners_.emplace(key, std::move(listener));
    set_stream_active(process, type, true);
  }

  void remove_io_listener(struct ged_subprocess* process, bu_process_io_t type)
  {
    const IoKey key{process, static_cast<int>(type)};
    auto entry = io_listeners_.find(key);
    if (entry == io_listeners_.end()) {
      set_stream_active(process, type, false);
      return;
    }

    stop_io_listener(*entry->second);
    io_listeners_.erase(entry);
    set_stream_active(process, type, false);
  }

  static void stop_io_listener(IoListener& listener)
  {
    if (!listener.watching) {
      return;
    }
    if (listener.polling) {
      Fl::remove_timeout(io_poll, &listener);
    } else {
      Fl::remove_fd(listener.descriptor, FL_READ);
    }
    listener.watching = false;
  }

  static bool process_polling_enabled()
  {
#if defined(_WIN32)
    return true;
#else
    const char* requested = std::getenv("FLTKGED_FORCE_PROCESS_POLL");
    return requested && !BU_STR_EQUAL(requested, "0");
#endif
  }

  void invoke_io_callback(struct ged* callback_ged, ged_io_func_t callback,
                          void* callback_data, int type)
  {
    if (!callback_ged || !callback) {
      return;
    }
    const std::size_t before = bu_vls_strlen(callback_ged->ged_result_str);
    try {
      callback(callback_data, type);
    } catch (const std::exception& error) {
      report_error(std::string("GED subprocess callback failed: ") +
                   error.what());
      return;
    } catch (...) {
      report_error("GED subprocess callback failed");
      return;
    }

    const std::size_t after = bu_vls_strlen(callback_ged->ged_result_str);
    if (after > before) {
      const std::size_t available = after - before;
      const std::size_t accepted =
        std::min(available, terminal_->limits().max_append_bytes);
      terminal_->append(std::string(bu_vls_cstr(callback_ged->ged_result_str) +
                                      before,
                                    accepted));
      if (accepted != available) {
        terminal_->append(kOutputTruncation);
        report_error("GED subprocess output exceeded the terminal limit");
      }
      if (before <= static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        bu_vls_trunc(callback_ged->ged_result_str, static_cast<int>(before));
      } else {
        bu_vls_trunc(callback_ged->ged_result_str, 0);
        report_error("GED subprocess result exceeded the supported size");
      }
    }
    views_->invalidate();
  }

  static void io_ready(int, void* data)
  {
    auto* listener = static_cast<IoListener*>(data);
    Impl* owner = listener->owner;
    struct ged* callback_ged = listener->process->gedp;
    const auto callback = listener->callback;
    void* callback_data = listener->callback_data;
    const int type = static_cast<int>(listener->type);
    owner->invoke_io_callback(callback_ged, callback, callback_data, type);
  }

  static void io_poll(void* data)
  {
    auto* listener = static_cast<IoListener*>(data);
    if (!listener || !listener->owner || !listener->process ||
        !listener->process->p) {
      return;
    }
    Impl* owner = listener->owner;
    const IoKey key{listener->process, static_cast<int>(listener->type)};
    const bool ready = bu_process_pending(listener->descriptor) != 0;
    const bool complete = bu_process_poll(listener->process->p, nullptr) != 0;
    if (ready || complete) {
      struct ged* callback_ged = listener->process->gedp;
      const auto callback = listener->callback;
      void* callback_data = listener->callback_data;
      const int type = static_cast<int>(listener->type);
      owner->invoke_io_callback(callback_ged, callback, callback_data, type);
    }

    const auto entry = owner->io_listeners_.find(key);
    if (entry != owner->io_listeners_.end() &&
        entry->second.get() == listener && listener->watching) {
      Fl::repeat_timeout(kProcessPollIntervalSeconds, io_poll, listener);
    }
  }

  void clear_io_listeners()
  {
    for (const auto& entry : io_listeners_) {
      stop_io_listener(*entry.second);
    }
    io_listeners_.clear();
  }

  static void set_stream_active(struct ged_subprocess* process,
                                bu_process_io_t type, bool active)
  {
    if (!process) {
      return;
    }
    const int value = active ? 1 : 0;
    switch (type) {
      case BU_PROCESS_STDIN: process->stdin_active = value; break;
      case BU_PROCESS_STDOUT: process->stdout_active = value; break;
      case BU_PROCESS_STDERR: process->stderr_active = value; break;
    }
  }

  static void open_fb_client(struct fbserv_obj* server, int client_index,
                             void*)
  {
    if (!server || !server->fbs_clientData || client_index < 0 ||
        client_index >= MAX_CLIENTS) {
      return;
    }
    auto& self = *static_cast<Impl*>(server->fbs_clientData);
    auto* client = &server->fbs_clients[client_index];
    if (client->fbsc_fd < 0) {
      self.report_error("The framebuffer IPC channel has no readable descriptor");
      return;
    }

    self.remove_fb_watch(client);
    auto watch = std::make_unique<FbWatch>();
    watch->owner = &self;
    watch->server = server;
    watch->client_index = client_index;
    watch->descriptor = client->fbsc_fd;
    watch->watching = true;
    auto* raw = watch.get();
    client->fbsc_chan = raw;
    self.fb_watches_.emplace(client, std::move(watch));
    Fl::add_fd(raw->descriptor, FL_READ, fb_client_ready, raw);
  }

  static void close_fb_client(struct fbserv_obj* server, int client_index)
  {
    if (!server || client_index < 0 || client_index >= MAX_CLIENTS) {
      return;
    }
    auto* client = &server->fbs_clients[client_index];
    auto* watch = static_cast<FbWatch*>(client->fbsc_chan);
    client->fbsc_chan = nullptr;
    if (watch && watch->owner) {
      watch->owner->remove_fb_watch(client);
    }
  }

  static void close_fb_server(struct fbserv_obj*)
  {
  }

  static void framebuffer_changed(void* data)
  {
    auto* self = static_cast<Impl*>(data);
    if (self && self->views_) {
      self->views_->invalidate();
    }
  }

  static void fb_client_ready(int, void* data)
  {
    auto* watch = static_cast<FbWatch*>(data);
    if (!watch || !watch->owner || !watch->server ||
        watch->client_index < 0 || watch->client_index >= MAX_CLIENTS) {
      return;
    }
    Impl* owner = watch->owner;
    struct fbserv_obj* server = watch->server;
    auto* client = &server->fbs_clients[watch->client_index];
    const int descriptor = watch->descriptor;
    if (watch->watching) {
      Fl::remove_fd(descriptor, FL_READ);
      watch->watching = false;
    }

    fbs_existing_client_handler(client, 0);
    if (owner->views_) {
      owner->views_->invalidate();
    }

    auto entry = owner->fb_watches_.find(client);
    if (entry != owner->fb_watches_.end() &&
        entry->second.get() == watch && client->fbsc_pkg &&
        client->fbsc_chan == watch) {
      watch->watching = true;
      Fl::add_fd(descriptor, FL_READ, fb_client_ready, watch);
    }
  }

  void remove_fb_watch(struct fbserv_client* client)
  {
    auto entry = fb_watches_.find(client);
    if (entry == fb_watches_.end()) {
      return;
    }
    if (entry->second->watching) {
      Fl::remove_fd(entry->second->descriptor, FL_READ);
    }
    fb_watches_.erase(entry);
  }

  void close_framebuffer_server()
  {
    if (ged_ && ged_->ged_fbs &&
        (ged_->ged_fbs->fbs_listener.fbsl_ipc_child || !fb_watches_.empty())) {
      fbs_close(ged_->ged_fbs);
    }
    for (const auto& entry : fb_watches_) {
      if (entry.second->watching) {
        Fl::remove_fd(entry.second->descriptor, FL_READ);
      }
    }
    fb_watches_.clear();
  }

  static void open_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    Fl_Native_File_Chooser chooser;
    chooser.title("Open BRL-CAD database");
    chooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
    chooser.filter("BRL-CAD database\t*.g");
    if (chooser.show() == 0 && chooser.filename()) {
      self.open_database(chooser.filename());
    }
  }

  static void new_database_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    Fl_Native_File_Chooser chooser;
    chooser.title("Create BRL-CAD database");
    chooser.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
    chooser.filter("BRL-CAD database\t*.g");
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }
    const std::filesystem::path path(chooser.filename());
    std::error_code error;
    if (std::filesystem::exists(path, error)) {
      self.report_panel_error("Refusing to overwrite an existing database: " +
                              path.string());
      fl_alert("The selected file already exists. Choose a new database name.");
      return;
    }
    if (error) {
      self.report_panel_error("Unable to inspect the selected database path: " +
                              error.message());
      return;
    }
    struct db_i* database = db_create(path.string().c_str(), 5);
    if (!database) {
      self.report_panel_error("Unable to create BRL-CAD database: " +
                              path.string());
      return;
    }
    db_close(database);
    self.open_database(path.string());
  }

  static void import_database_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->dbip) {
      self.report_panel_error("Open a destination database before importing");
      return;
    }
    Fl_Native_File_Chooser chooser;
    chooser.title("Import BRL-CAD database");
    chooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
    chooser.filter("BRL-CAD database\t*.g");
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }
    const char* response = fl_input(
      "Prefix for imported names (use / for automatic conflict renaming):",
      "/");
    if (!response) {
      return;
    }
    const std::string affix = *response ? response : "/";
    self.execute("concat -p " + quote(chooser.filename()) + " " +
                 quote(affix));
  }

  static void export_selected_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.selected_object_.empty()) {
      self.report_panel_error("Select an object to export first");
      return;
    }
    Fl_Native_File_Chooser chooser;
    chooser.title("Export selected object");
    chooser.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
    chooser.filter("BRL-CAD database\t*.g");
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }
    const std::filesystem::path path(chooser.filename());
    std::error_code error;
    if (std::filesystem::exists(path, error)) {
      self.report_panel_error("Refusing to append to an existing export: " +
                              path.string());
      fl_alert("The selected file already exists. Choose a new export name.");
      return;
    }
    if (error) {
      self.report_panel_error("Unable to inspect the export path: " +
                              error.message());
      return;
    }
    self.execute("keep " + quote(path.string()) + " " +
                 quote(self.selected_object_));
  }

  static void source_tcl_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.tcl_->available()) {
      self.report_panel_error(self.tcl_->last_error());
      fl_alert("%s", self.tcl_->last_error().c_str());
      return;
    }
    Fl_Native_File_Chooser chooser;
    chooser.title("Load Tcl script");
    chooser.type(Fl_Native_File_Chooser::BROWSE_FILE);
    chooser.filter("Tcl script\t*.tcl\nAll files\t*");
    if (chooser.show() != 0 || !chooser.filename()) {
      return;
    }
    self.activate_view();
    const int status = self.tcl_->source(chooser.filename());
    const std::string& output = self.tcl_->result();
    if (!output.empty()) {
      self.terminal_->append(output);
      if (output.back() != '\n') {
        self.terminal_->append("\n");
      }
    }
    if (status != 0) {
      self.report_panel_error(self.tcl_->last_error());
    }
  }

  static void enable_tk_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.tcl_->enable_legacy_tk()) {
      self.report_panel_error(self.tcl_->last_error());
      fl_alert("%s", self.tcl_->last_error().c_str());
      return;
    }
    self.terminal_->append(
      "Legacy Tk event integration enabled (Tk " + self.tcl_->result() +
      "). Native FLTK panels remain available through the fltkged Tcl command.\n");
  }

  static void quit_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->request_quit();
  }

  static void redraw_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->views_->invalidate();
  }

  static void autoview_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("autoview");
  }

  static void refresh_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->refresh_database_panels();
  }

  static void preferences_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.preferences_editor_) {
      self.preferences_editor_ = std::make_unique<Fl_Ged_Preferences>();
      self.preferences_editor_->changed_handler(
        [&self](const Fl_Ged_Preferences::Settings& settings) {
          self.apply_preferences(settings, true);
        });
    }
    self.preferences_editor_->settings(self.preferences_settings_);
    self.preferences_editor_->show();
  }

  static void edit_primitive_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.selected_object_.empty()) {
      self.report_panel_error("Select a primitive in the hierarchy first");
      fl_alert("Select a primitive in the hierarchy first");
      return;
    }
    if (!self.primitive_editor_) {
      self.primitive_editor_ = std::make_unique<Fl_Ged_Primitive_Editor>(
        self.views_->active_view().preference());
      self.primitive_editor_->changed_handler([&self]() {
        self.refresh_database_panels();
        self.views_->invalidate();
      });
      self.primitive_editor_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    if (!self.primitive_editor_->open(self.ged_, self.selected_object_)) {
      fl_alert("%s", self.primitive_editor_->last_error().c_str());
      return;
    }
    self.primitive_editor_->show();
  }

  static void create_primitive_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->dbip) {
      self.report_panel_error("Open a database before creating a primitive");
      fl_alert("Open a database before creating a primitive");
      return;
    }
    if (!self.create_primitive_) {
      self.create_primitive_ = std::make_unique<Fl_Ged_Create_Primitive>();
      self.create_primitive_->created_handler(
        [&self](const std::string& object_name) {
          self.selected_object_ = object_name;
          self.refresh_database_panels();
          self.attributes_->object(object_name);
          self.views_->invalidate();
          self.refresh_faceplate();
        });
      self.create_primitive_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    self.create_primitive_->ged(self.ged_);
    self.create_primitive_->show();
  }

  static void search_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->dbip) {
      self.report_panel_error("Open a database before searching");
      fl_alert("Open a database before searching");
      return;
    }
    if (!self.search_) {
      self.search_ = std::make_unique<Fl_Ged_Search>();
      self.search_->selection_handler([&self](const std::string& object_name) {
        self.selected_object_ = object_name;
        self.attributes_->object(object_name);
        self.refresh_faceplate();
      });
      self.search_->command_handler(
        [&self](const std::string& command) { self.execute(command); });
      self.search_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    self.search_->ged(self.ged_);
    self.search_->show();
  }

  static void edit_combination_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.selected_object_.empty()) {
      self.report_panel_error("Select a combination in the hierarchy first");
      fl_alert("Select a combination in the hierarchy first");
      return;
    }
    if (!self.combination_editor_) {
      self.combination_editor_ =
        std::make_unique<Fl_Ged_Combination_Editor>();
      self.combination_editor_->changed_handler([&self]() {
        self.refresh_database_panels();
        self.views_->invalidate();
      });
      self.combination_editor_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    if (!self.combination_editor_->open(self.ged_, self.selected_object_)) {
      fl_alert("%s", self.combination_editor_->last_error().c_str());
      return;
    }
    self.combination_editor_->show();
  }

  static void edit_sketch_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.selected_object_.empty()) {
      self.report_panel_error("Select a sketch in the hierarchy first");
      fl_alert("Select a sketch in the hierarchy first");
      return;
    }
    if (!self.sketch_editor_) {
      self.sketch_editor_ = std::make_unique<Fl_Ged_Sketch_Editor>(
        self.views_->active_view().preference());
      self.sketch_editor_->changed_handler([&self]() {
        self.refresh_database_panels();
        self.views_->invalidate();
      });
      self.sketch_editor_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    if (!self.sketch_editor_->open(self.ged_, self.selected_object_)) {
      fl_alert("%s", self.sketch_editor_->last_error().c_str());
      return;
    }
    self.sketch_editor_->show();
  }

  static void transform_editor_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (self.selected_path_.empty()) {
      self.report_panel_error("Select an object in the hierarchy or display first");
      return;
    }
    if (!self.transform_editor_) {
      self.transform_editor_ = std::make_unique<Fl_Ged_Transform_Editor>();
      self.transform_editor_->command_handler(
        [&self](const std::string& command) { return self.execute(command); });
      self.transform_editor_->changed_handler([&self]() {
        self.refresh_database_panels();
        self.views_->invalidate();
      });
      self.transform_editor_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    self.transform_editor_->ged(self.ged_);
    self.transform_editor_->selection(self.selected_path_);
    self.transform_editor_->show();
  }

  static void single_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->views_->layout(Fl_Cad_Quad_View::Layout::Single);
  }

  static void render_tools_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.render_tools_) {
      self.render_tools_ = std::make_unique<Fl_Ged_Render_Tools>();
      self.render_tools_->command_handler(
        [&self](const std::string& command) { return self.execute(command); });
      self.render_tools_->error_handler([&self](const std::string& message) {
        self.report_panel_error(message);
      });
    }
    self.render_tools_->show();
  }

  static void center_ray_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("vquery_ray 0 0");
  }

  static void abort_raytrace_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("rtabort");
  }

  static bool require_selection(Impl& self)
  {
    if (!self.selected_object_.empty()) {
      return true;
    }
    self.report_panel_error("Select a database object first");
    return false;
  }

  static void analyze_selected_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (require_selection(self)) {
      self.execute("analyze " + quote(self.selected_object_));
    }
  }

  static void gqa_selected_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (require_selection(self)) {
      self.execute("gqa " + quote(self.selected_object_));
    }
  }

  static void check_overlaps_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (require_selection(self)) {
      self.execute("check overlaps -i -o " + quote(self.selected_object_));
    }
  }

  static void specialist_tools_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->dbip) {
      self.report_panel_error("Open a database before using specialist tools");
      return;
    }
    if (!self.specialist_tools_) {
      self.specialist_tools_ = std::make_unique<Fl_Ged_Specialist_Tools>();
      self.specialist_tools_->command_handler(
        [&self](const std::string& command) { return self.execute(command); });
      self.specialist_tools_->changed_handler([&self]() {
        self.refresh_database_panels();
        self.views_->invalidate();
      });
      self.specialist_tools_->error_handler(
        [&self](const std::string& message) { self.report_panel_error(message); });
    }
    self.specialist_tools_->ged(self.ged_);
    self.specialist_tools_->selection(self.selected_object_);
    self.specialist_tools_->show();
  }

  static void quad_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->views_->layout(Fl_Cad_Quad_View::Layout::Quad);
  }

  static void top_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("ae 270 90");
  }

  static void front_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("ae 0 0");
  }

  static void right_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("ae 270 0");
  }

  static void left_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("ae 90 0");
  }

  static void aet_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->execute("ae 35 25");
  }

  static void add_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->add_stored_view();
  }

  static void next_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->step_stored_view(1);
  }

  static void previous_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->step_stored_view(-1);
  }

  static void delete_view_menu(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->delete_stored_view();
  }

  static void documentation_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.help_) {
      self.help_ = std::make_unique<Fl_Ged_Help>();
    }
    if (!self.help_->load_default()) {
      self.report_panel_error(self.help_->last_error());
      fl_alert("%s", self.help_->last_error().c_str());
      return;
    }
    self.help_->show();
  }

  static void about_menu(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.about_ = std::make_unique<Fl_Ged_About>(
      backend_name(self.views_->active_view().backend()));
    self.about_->show();
  }

  FltkGedWindow& owner_;
  Fl_Menu_Bar* menu_ = nullptr;
  Fl_Tile* main_tile_ = nullptr;
  Fl_Tabs* panels_ = nullptr;
  Fl_Group* hierarchy_group_ = nullptr;
  Fl_Group* attributes_group_ = nullptr;
  Fl_Group* view_group_ = nullptr;
  Fl_Group* interaction_group_ = nullptr;
  Fl_Ged_Hierarchy* hierarchy_ = nullptr;
  Fl_Ged_Attributes* attributes_ = nullptr;
  Fl_Ged_View_Controls* view_controls_ = nullptr;
  Fl_Ged_Interaction* interaction_ = nullptr;
  std::unique_ptr<Fl_Double_Window> graphics_window_;
  std::unique_ptr<Fl_Double_Window> panels_window_;
  Fl_Menu_Bar* graphics_menu_ = nullptr;
  bool combined_ = false;
  Fl_Cad_Quad_View* views_ = nullptr;
  Fl_Cad_Terminal* terminal_ = nullptr;
  struct ged* ged_ = nullptr;
  struct db_i* database_before_open_ = nullptr;
  std::unique_ptr<Fl_Ged_Tcl> tcl_;
  FltkGedWindow::CommandLanguage command_language_ = FltkGedWindow::CommandLanguage::Tcl;
  std::string pending_script_;
  bool pending_tcl_more_ = false;
  bool terminal_submission_ = false;
  std::vector<std::string> pending_arguments_;
  std::vector<bool> draw_autoview_;
  std::optional<std::size_t> pending_history_begin_;
  using IoKey = std::pair<struct ged_subprocess*, int>;
  std::map<IoKey, std::unique_ptr<IoListener>> io_listeners_;
  std::map<struct fbserv_client*, std::unique_ptr<FbWatch>> fb_watches_;
  std::array<Faceplate, 4> faceplates_;
  std::array<std::vector<ViewSnapshot>, 4> stored_views_;
  std::array<std::size_t, 4> stored_view_indices_{{
    std::numeric_limits<std::size_t>::max(),
    std::numeric_limits<std::size_t>::max(),
    std::numeric_limits<std::size_t>::max(),
    std::numeric_limits<std::size_t>::max()}};
  std::unique_ptr<Fl_Ged_Help> help_;
  std::unique_ptr<Fl_Ged_About> about_;
  std::unique_ptr<Fl_Ged_Primitive_Editor> primitive_editor_;
  std::unique_ptr<Fl_Ged_Preferences> preferences_editor_;
  std::unique_ptr<Fl_Ged_Combination_Editor> combination_editor_;
  std::unique_ptr<Fl_Ged_Sketch_Editor> sketch_editor_;
  std::unique_ptr<Fl_Ged_Create_Primitive> create_primitive_;
  std::unique_ptr<Fl_Ged_Search> search_;
  std::unique_ptr<Fl_Ged_Render_Tools> render_tools_;
  std::unique_ptr<Fl_Ged_Transform_Editor> transform_editor_;
  std::unique_ptr<Fl_Ged_Specialist_Tools> specialist_tools_;
  std::mutex log_mutex_;
  std::string pending_log_;
  bool log_truncated_ = false;
  std::string selected_object_;
  std::string selected_path_;
  std::string last_error_;
  Fl_Ged_Preferences::Settings preferences_settings_;
};

FltkGedWindow::FltkGedWindow(int width, int height,
                             Fl_Cad_View::Backend backend, const char* title)
  : Fl_Double_Window(width, height, title),
    impl_(std::make_unique<Impl>(*this, backend))
{
}

FltkGedWindow::~FltkGedWindow() { hide(); }

bool FltkGedWindow::open_database(const std::string& path)
{
  return impl_->open_database(path);
}

int FltkGedWindow::execute(const std::string& command)
{
  return impl_->execute(command);
}
int FltkGedWindow::execute_prompt(const std::string& command)
{
  return impl_->execute_prompt(command);
}
void FltkGedWindow::command_language(CommandLanguage language)
{
  impl_->command_language(language);
}
FltkGedWindow::CommandLanguage FltkGedWindow::command_language() const noexcept
{
  return impl_->command_language_;
}

struct ged* FltkGedWindow::gedp() noexcept { return impl_->ged_; }
Fl_Cad_View& FltkGedWindow::view() noexcept
{
  return impl_->views_->active_view();
}
Fl_Cad_Terminal& FltkGedWindow::terminal() noexcept { return *impl_->terminal_; }
const std::string& FltkGedWindow::last_error() const noexcept
{
  return impl_->last_error_;
}

void FltkGedWindow::show()
{
  Fl_Double_Window::show();
  if (!impl_->combined_) impl_->graphics_window_->show();
  impl_->terminal_->output().take_focus();
}

void FltkGedWindow::hide()
{
  if (impl_) impl_->hide_windows();
  Fl_Double_Window::hide();
}

void FltkGedWindow::combined(bool value) { impl_->combined(value); }
bool FltkGedWindow::combined() const noexcept { return impl_->combined_; }
Fl_Window& FltkGedWindow::graphics_window() noexcept { return *impl_->graphics_window_; }

void FltkGedWindow::resize(int x, int y, int width, int height)
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
