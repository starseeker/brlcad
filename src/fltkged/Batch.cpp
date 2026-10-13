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

#include "Batch.H"
#include "fltkcad/Fl_Ged_Tcl.H"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

extern "C" {
#include "bu/str.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/commands.h"
#include "dm.h"
}

namespace {
constexpr std::size_t kMaxCommandBytes = 4U * 1024U * 1024U;
constexpr std::size_t kMaxArguments = 4096;

struct GedCloser {
  void operator()(struct ged* context) const
  {
    if (!context) return;
    ged_subprocesses_terminate(context);
    const char* arguments[] = {"Z"};
    ged_exec_zap(context, 1, arguments);
    if (context->ged_gvp && context->ged_gvp->dmp) {
      dm_close(static_cast<struct dm*>(context->ged_gvp->dmp));
      context->ged_gvp->dmp = nullptr;
    }
    ged_close(context);
  }
};

void activate_scene_object(void*, struct bv_scene_obj* object)
{
  if (object) object->s_flag = UP;
}
}

int fltkged_batch(const std::string& database, const std::string& command,
                   bool raw_ged)
{
  std::unique_ptr<struct ged, GedCloser> context(
    database.empty() ? ged_create() : ged_open("db", database.c_str(), 0));
  if (!context) {
    std::cerr << "Unable to open database: " << database << '\n';
    return 1;
  }
  context->ged_create_vlist_scene_obj_callback = activate_scene_object;
  context->ged_gvp->dmp = dm_open(nullptr, nullptr, "nu", 0, nullptr);
  if (!context->ged_gvp->dmp) {
    std::cerr << "Unable to initialize the batch display manager\n";
    return 1;
  }
  const vect_t orientation = {35.0, 25.0, 0.0};
  bv_view_set_aet(context->ged_gvp, orientation);
  Fl_Ged_Tcl interpreter;
  interpreter.ged(context.get());
  if (!raw_ged && !interpreter.available()) {
    std::cerr << interpreter.last_error() << '\n';
    return 1;
  }

  bool quit = false;
  auto evaluate = [&](const std::string& script) {
    if (script.size() > kMaxCommandBytes || script.find('\0') != std::string::npos) {
      std::cerr << "Command exceeds the configured input limit\n";
      return 1;
    }
    if (!raw_ged) {
      const int status = interpreter.evaluate(script);
      quit = interpreter.exit_requested();
      if (quit) return interpreter.exit_status();
      const std::string& output = status == 0 ? interpreter.result() : interpreter.last_error();
      if (!output.empty()) {
        (status == 0 ? std::cout : std::cerr) << output;
        if (output.back() != '\n') (status == 0 ? std::cout : std::cerr) << '\n';
      }
      return status == 0 ? 0 : 1;
    }
    std::vector<char> text(script.begin(), script.end());
    text.push_back('\0');
    std::vector<char*> arguments(kMaxArguments + 1, nullptr);
    const int count = bu_argv_from_string(arguments.data(), static_cast<int>(arguments.size()), text.data());
    if (!count) return 0;
    if (static_cast<std::size_t>(count) > kMaxArguments) {
      std::cerr << "GED argument limit exceeded\n";
      return 1;
    }
    const int status = ged_edit_session_exec(context.get(), count, const_cast<const char**>(arguments.data()));
    const std::string output = bu_vls_cstr(context->ged_result_str);
    bu_vls_trunc(context->ged_result_str, 0);
    if (!output.empty()) std::cout << output << (output.back() == '\n' ? "" : "\n");
    quit = (status & GED_EXIT) != 0;
    return (status & ~GED_EXIT) == BRLCAD_OK ? 0 : 1;
  };

  if (!command.empty()) return evaluate(command);
  std::string script;
  std::string line;
  int result = 0;
  while (!quit && std::getline(std::cin, line)) {
    if (script.size() + line.size() + 1 > kMaxCommandBytes) {
      std::cerr << "Command exceeds the configured input limit\n";
      return 1;
    }
    script += line;
    script += '\n';
    if (!raw_ged && !interpreter.command_complete(script)) continue;
    const int status = evaluate(script);
    if (status) result = status;
    script.clear();
  }
  if (!script.empty()) {
    std::cerr << "Incomplete Tcl command at end of input\n";
    return 1;
  }
  return quit ? interpreter.exit_status() : result;
}
