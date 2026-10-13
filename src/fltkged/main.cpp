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

#include "FltkGedWindow.H"
#include "Batch.H"
#include "Fl_Ged_Localization.H"

#include "fltkcad/Fl_Cad_Terminal.H"

#include <FL/Fl.H>
#include <FL/Fl_Preferences.H>
#include <FL/Fl_Terminal.H>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>

extern "C" {
#include "bu/app.h"
#include "ged.h"
}

namespace {

bool has_visible_rgb(const Fl_Cad_View::Image& image)
{
  for (std::size_t offset = 0; offset + 3 < image.rgba.size(); offset += 4) {
    if (image.rgba[offset] != 0 || image.rgba[offset + 1] != 0 ||
        image.rgba[offset + 2] != 0) {
      return true;
    }
  }
  return false;
}

bool images_differ(const Fl_Cad_View::Image& first,
                   const Fl_Cad_View::Image& second)
{
  return first.width != second.width || first.height != second.height ||
         first.rgba != second.rgba;
}

void usage(const char* program)
{
  std::cerr << "Usage: " << program
            << " [-c] [--shell=tcl|ged] [--software|--opengl] [--combined]"
               " [--smoke|--render-smoke] [database.g [command ...]]\n";
}

} // namespace

int main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  fltkged::initialize_locale();
  Fl_Cad_View::Backend backend = Fl_Cad_View::Backend::Auto;
  bool combined = false;
  bool smoke = false;
  bool render_smoke = false;
  bool batch = false;
  bool raw_ged = false;
  bool require_tcl = false;
  std::string database;
  std::string command;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "-c") {
      batch = true;
    } else if (argument == "--shell=ged") {
      raw_ged = true;
      require_tcl = false;
    } else if (argument == "--shell=tcl") {
      raw_ged = false;
      require_tcl = true;
    } else if (argument == "--help" || argument == "-h") {
      usage(argv[0]);
      return EXIT_SUCCESS;
    } else if (argument == "--software") {
      backend = Fl_Cad_View::Backend::OSMesa;
    } else if (argument == "--opengl") {
      backend = Fl_Cad_View::Backend::SystemOpenGL;
    } else if (argument == "--combined") {
      combined = true;
    } else if (argument == "--smoke") {
      smoke = true;
    } else if (argument == "--render-smoke") {
      smoke = true;
      render_smoke = true;
    } else if (!argument.empty() && argument.front() == '-') {
      usage(argv[0]);
      return EXIT_FAILURE;
    } else if (database.empty()) {
      database = argument;
      for (++index; index < argc; ++index) {
        if (!command.empty()) command += ' ';
        command += argv[index];
      }
      break;
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (batch) return fltkged_batch(database, command, raw_ged);

  if (smoke && std::getenv("DISPLAY") == nullptr) {
    return 77;
  }
  if (smoke) {
    Fl_Preferences::file_access(Fl_Preferences::NONE);
  }

  FltkGedWindow window(900, combined ? 900 : 340, backend);
  if (require_tcl && window.command_language() != FltkGedWindow::CommandLanguage::Tcl) {
    std::cerr << "Tcl is unavailable; use --shell=ged for this build\n";
    return EXIT_FAILURE;
  }
  if (raw_ged) window.command_language(FltkGedWindow::CommandLanguage::Ged);
  window.combined(combined);
  if (!database.empty() && !window.open_database(database)) {
    std::cerr << window.last_error() << '\n';
    return EXIT_FAILURE;
  }
  window.show();
  if (!command.empty() && window.execute_prompt(command) != BRLCAD_OK)
    return EXIT_FAILURE;

  if (!smoke) {
    return Fl::run();
  }

  const auto initialization_deadline = std::chrono::steady_clock::now() +
                                       std::chrono::seconds(2);
  while (!window.view().initialized() &&
         std::chrono::steady_clock::now() < initialization_deadline) {
    Fl::wait(0.01);
  }
  if (!window.view().initialize()) {
    std::cerr << window.view().last_error() << '\n';
    return EXIT_FAILURE;
  }
  if (!database.empty()) {
    if (window.execute("draw all.g") != 0 || window.execute("autoview") != 0) {
      std::cerr << "GED smoke command failed\n";
      return EXIT_FAILURE;
    }
  } else if (window.execute("help") != 0) {
    std::cerr << "GED help command failed\n";
    return EXIT_FAILURE;
  }
  window.view().invalidate();
  const auto draw_deadline = std::chrono::steady_clock::now() +
                             std::chrono::milliseconds(250);
  do {
    Fl::wait(0.01);
  } while (std::chrono::steady_clock::now() < draw_deadline);
  auto capture = window.view().capture();
  if (!capture || capture->rgba.size() !=
                    static_cast<std::size_t>(capture->width) *
                      static_cast<std::size_t>(capture->height) * 4U) {
    std::cerr << window.view().last_error() << '\n';
    return EXIT_FAILURE;
  }
  if (!database.empty() && !has_visible_rgb(*capture)) {
    std::cerr << "Rendered database is blank\n";
    return EXIT_FAILURE;
  }
  if (render_smoke) {
    if (database.empty()) {
      std::cerr << "Render smoke testing requires a database\n";
      return EXIT_FAILURE;
    }
    const Fl_Cad_View::Image before = *capture;
    if (window.execute("ert -P 1") != 0) {
      std::cerr << "Unable to start embedded raytrace: "
                << window.last_error() << '\n';
      if (const char* output = window.terminal().output().text(); output) {
        std::cerr << output << '\n';
      }
      return EXIT_FAILURE;
    }
    bool updated = false;
    const auto render_deadline = std::chrono::steady_clock::now() +
                                 std::chrono::seconds(12);
    while (!updated && std::chrono::steady_clock::now() < render_deadline) {
      Fl::wait(0.02);
      auto rendered = window.view().capture();
      updated = rendered && images_differ(before, *rendered);
    }
    if (!updated) {
      window.execute("rtabort");
      std::cerr << "Embedded raytrace produced no framebuffer update\n";
      if (const char* output = window.terminal().output().text(); output) {
        std::cerr << output << '\n';
      }
      return EXIT_FAILURE;
    }
  }
  window.hide();
  return EXIT_SUCCESS;
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
