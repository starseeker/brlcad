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

#include "fltkcad/Fl_Ged_Tcl.H"

#include <FL/Fl.H>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bu/app.h"
#include "bu/process.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/commands.h"
}

#if defined(FLTKCAD_HAVE_TCL)
#  include <tcl.h>
#  include "tclcad/setup.h"
#endif

namespace {

constexpr std::size_t kMaxScriptBytes = 4U * 1024U * 1024U;
constexpr std::size_t kMaxResultBytes = 4U * 1024U * 1024U;
constexpr int kMaxGedArguments = 4096;
constexpr int kEventPollMilliseconds = 10;
constexpr double kTclEventIntervalSeconds = kEventPollMilliseconds / 1000.0;
constexpr int kMaxTclEventsPerTick = 32;

} // namespace

class Fl_Ged_Tcl::Impl
{
public:
#if defined(FLTKCAD_HAVE_TCL)
  struct Command
  {
    Impl* owner = nullptr;
    std::string name;
  };
#endif

  Impl()
  {
#if defined(FLTKCAD_HAVE_TCL)
    Tcl_FindExecutable(bu_getprogname());
    interpreter_ = Tcl_CreateInterp();
    if (!interpreter_) {
      last_error_ = "Unable to create a Tcl interpreter";
      return;
    }
    if (tclcad_init(interpreter_, 0, nullptr) != TCL_OK ||
        Tcl_EvalEx(interpreter_, "auto_load help; set ::glob_compat_mode 0",
                   -1, TCL_EVAL_GLOBAL) != TCL_OK) {
      last_error_ = Tcl_GetStringResult(interpreter_);
      Tcl_DeleteInterp(interpreter_);
      interpreter_ = nullptr;
      return;
    }
    Tcl_CreateObjCommand(interpreter_, "ged", ged_dispatch, this, nullptr);
    Tcl_CreateObjCommand(interpreter_, "fltkged", ui_dispatch, this, nullptr);
    Tcl_EvalEx(interpreter_, "namespace eval ::ged {}", -1, TCL_EVAL_GLOBAL);
    Tcl_EvalEx(interpreter_,
               "namespace eval ::fltkged {}; package provide fltkged 1.0",
               -1, TCL_EVAL_GLOBAL);
    Tcl_EvalEx(interpreter_,
               "rename ::puts ::fltkged::tcl_puts; rename ::exit ::fltkged::tcl_exit",
               -1, TCL_EVAL_GLOBAL);
    Tcl_CreateObjCommand(interpreter_, "puts", puts_dispatch, this, nullptr);
    Tcl_CreateObjCommand(interpreter_, "exit", exit_dispatch, this, nullptr);

    const char* const* names = nullptr;
    const std::size_t count = ged_cmd_list(&names);
    commands_.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
      if (!names[index] || !*names[index]) {
        continue;
      }
      auto command = std::make_unique<Command>();
      command->owner = this;
      command->name = names[index];
      Command* data = command.get();
      Tcl_CreateObjCommand(interpreter_, ("::ged::" + command->name).c_str(),
                           direct_dispatch, data, nullptr);
      Tcl_CreateObjCommand(interpreter_, ("::_mged_" + command->name).c_str(),
                           direct_dispatch, data, nullptr);
      Tcl_CmdInfo existing;
      if (!Tcl_GetCommandInfo(interpreter_, command->name.c_str(), &existing)) {
        Tcl_CreateObjCommand(interpreter_, command->name.c_str(),
                             direct_dispatch, data, nullptr);
      }
      commands_.push_back(std::move(command));
    }
    Tcl_CreateObjCommand(interpreter_, "::ged::commands", commands_dispatch, this, nullptr);
    Tcl_CreateObjCommand(interpreter_, "refresh", refresh_dispatch, this, nullptr);
    Tcl_CreateObjCommand(interpreter_, "output_hook", output_hook_dispatch, this, nullptr);
    Tcl_TraceVar(interpreter_, "perspective", TCL_GLOBAL_ONLY | TCL_TRACE_READS |
                 TCL_TRACE_WRITES, perspective_trace, this);
    for (const char* name : {"base2local", "local2base"})
      Tcl_TraceVar(interpreter_, name, TCL_GLOBAL_ONLY | TCL_TRACE_READS, units_trace, this);
    available_ = true;
#else
    last_error_ = "Tcl support was not available when libfltkcad was built";
#endif
  }

  ~Impl()
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (interpreter_) {
      bu_log_delete_hook(capture_output_hook, this);
      bind_ged(nullptr);
      event_processing(false);
      Tcl_DeleteInterp(interpreter_);
    }
#endif
  }

  int evaluate(const std::string& script)
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (!available_ || !interpreter_) {
      return fail(last_error_.empty() ? "Tcl is unavailable" : last_error_);
    }
    if (script.find('\0') != std::string::npos || script.size() > kMaxScriptBytes ||
        script.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      return fail("Tcl script exceeds the configured size limit");
    }
    if (evaluating_ || servicing_events_) {
      return fail("A Tcl command is already running");
    }
    evaluating_ = true;
    exit_requested_ = false;
    exit_status_ = 0;
    command_status_ = BRLCAD_OK;
    if (!set_limits()) {
      evaluating_ = false;
      return fail(last_error_);
    }
    Tcl_Obj* command = Tcl_NewStringObj(script.data(), static_cast<int>(script.size()));
    Tcl_IncrRefCount(command);
    const int status = Tcl_EvalObjEx(interpreter_, command, TCL_EVAL_GLOBAL);
    Tcl_DecrRefCount(command);
    clear_limits();
    evaluating_ = false;
    copy_result();
    start_tk_events_if_loaded();
    if (status != TCL_OK) {
      last_error_ = result_.empty() ? "Tcl evaluation failed" : result_;
    } else {
      last_error_.clear();
    }
    return status;
#else
    static_cast<void>(script);
    return fail(last_error_);
#endif
  }

  int source(const std::string& path)
  {
    if (path.empty() || path.find('\0') != std::string::npos) {
      return fail("A non-empty Tcl script path is required");
    }
    std::ifstream input(path, std::ios::binary);
    if (!input) {
      return fail("Unable to open Tcl script: " + path);
    }
    input.seekg(0, std::ios::end);
    const std::streamoff length = input.tellg();
    if (length < 0 || static_cast<unsigned long long>(length) > kMaxScriptBytes) {
      return fail("Tcl script exceeds the configured size limit");
    }
    input.seekg(0, std::ios::beg);
    std::string script(static_cast<std::size_t>(length), '\0');
    input.read(script.data(), length);
    if (input.gcount() != length) {
      return fail("Unable to read Tcl script: " + path);
    }
    char extra = '\0';
    if (input.get(extra)) {
      return fail("Tcl script changed while reading or exceeds the configured size limit");
    }
    return evaluate(script);
  }

  Fl_Ged_Tcl::ChangedHandler changed_handler_;
  Fl_Ged_Tcl::OutputHandler output_handler_;
  Fl_Ged_Tcl::UiHandler ui_handler_;
  struct ged* ged_ = nullptr;
  std::string result_;
  std::string last_error_;
  bool available_ = false;
  bool tk_active_ = false;
  bool events_enabled_ = false;
  bool evaluating_ = false;
  bool servicing_events_ = false;
  bool exit_requested_ = false;
  int exit_status_ = 0;
  int command_status_ = BRLCAD_OK;

  void interrupt()
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (interpreter_ && (evaluating_ || servicing_events_))
      Tcl_CancelEval(interpreter_, nullptr, nullptr, TCL_CANCEL_UNWIND);
#endif
  }

  void bind_ged(struct ged* value)
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (value == ged_) return;
    if (owns_process_io_) {
      ged_subprocesses_terminate(ged_);
      if (process_timer_) Tcl_DeleteTimerHandler(process_timer_);
      process_timer_ = nullptr;
      process_listeners_.clear();
      ged_->ged_create_io_handler = nullptr;
      ged_->ged_delete_io_handler = nullptr;
      ged_->ged_io_data = nullptr;
      owns_process_io_ = false;
    }
    ged_ = value;
    if (ged_ && !ged_->ged_create_io_handler && !ged_->ged_delete_io_handler) {
      ged_->ged_create_io_handler = create_process_io;
      ged_->ged_delete_io_handler = delete_process_io;
      ged_->ged_io_data = this;
      owns_process_io_ = true;
    }
#else
    ged_ = value;
#endif
  }

  bool command_complete(const std::string& script) const
  {
#if defined(FLTKCAD_HAVE_TCL)
    return Tcl_CommandComplete(script.c_str()) != 0;
#else
    static_cast<void>(script);
    return true;
#endif
  }

  void event_processing(bool enabled)
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (!interpreter_ || enabled == events_enabled_) return;
    events_enabled_ = enabled;
    if (enabled) {
      Tcl_CreateEventSource(setup_fltk_events, check_fltk_events, this);
      Fl::add_timeout(kTclEventIntervalSeconds, service_tcl_events, this);
    } else {
      Fl::remove_timeout(service_tcl_events, this);
      Tcl_DeleteEventSource(setup_fltk_events, check_fltk_events, this);
      Tcl_DeleteEvents(delete_fltk_event, this);
      fltk_event_queued_ = false;
    }
#else
    static_cast<void>(enabled);
#endif
  }

  bool enable_legacy_tk()
  {
#if defined(FLTKCAD_HAVE_TCL)
    if (!available_ || !interpreter_) {
      fail(last_error_.empty() ? "Tcl is unavailable" : last_error_);
      return false;
    }
    if (!set_limits()) return false;
    const char* version = Tcl_PkgRequire(interpreter_, "Tk", nullptr, 0);
    clear_limits();
    if (!version) {
      copy_result();
      last_error_ = result_.empty()
                      ? "Tk is not installed for the embedded Tcl interpreter"
                      : result_;
      return false;
    }
    start_tk_events_if_loaded();
    static constexpr const char* compatibility = R"tcl(
      if {![llength [info commands ::opendb]]} {
        proc ::opendb {args} {
          if {[llength $args]} {
            return -code error "opening databases from a legacy Tk script is not supported; use the FLTK File menu"
          }
          return [fltkged database]
        }
      }
      if {![llength [info commands ::_mged_opendb]]} {
        proc ::_mged_opendb {} { return [fltkged database] }
      }
      if {![llength [info commands ::winset]]} {
        proc ::winset {args} { return mged }
      }
      if {![llength [info commands ::cmd_win]]} {
        proc ::cmd_win {args} { return mged }
      }
      if {![llength [info commands ::mged_update]]} {
        proc ::mged_update {{non_blocking 0}} {
          update idletasks
          if {!$non_blocking} { update }
        }
      }
      set ::mged_players {mged}
      set ::mged_gui(mged,screen) [winfo screen .]
      set ::mged_gui(mged,active_dm) {}
      set ::mged_gui(mged,dmc) {}
      set ::mged_gui(mged,cmd) {}
    )tcl";
    if (Tcl_EvalEx(interpreter_, compatibility, -1, TCL_EVAL_GLOBAL) != TCL_OK) {
      copy_result();
      last_error_ = result_.empty()
                      ? "Unable to initialize the legacy MGED Tcl adapter"
                      : result_;
      return false;
    }
    result_ = version;
    last_error_.clear();
    return tk_active_;
#else
    fail(last_error_);
    return false;
#endif
  }

private:
  int fail(std::string message)
  {
    result_.clear();
    last_error_ = std::move(message);
#if defined(FLTKCAD_HAVE_TCL)
    return TCL_ERROR;
#else
    return 1;
#endif
  }

#if defined(FLTKCAD_HAVE_TCL)
  struct ProcessListener {
    struct ged_subprocess* process;
    bu_process_io_t stream;
    ged_io_func_t callback;
    void* data;
  };
  using ProcessKey = std::pair<struct ged_subprocess*, bu_process_io_t>;

  static void set_process_stream(struct ged_subprocess* process,
                                 bu_process_io_t stream, bool active)
  {
    switch (stream) {
      case BU_PROCESS_STDIN: process->stdin_active = active; break;
      case BU_PROCESS_STDOUT: process->stdout_active = active; break;
      case BU_PROCESS_STDERR: process->stderr_active = active; break;
    }
  }

  static void create_process_io(struct ged_subprocess* process, bu_process_io_t stream,
                                 ged_io_func_t callback, void* data)
  {
    auto& self = *static_cast<Impl*>(process->gedp->ged_io_data);
    self.process_listeners_[{process, stream}] = {process, stream, callback, data};
    set_process_stream(process, stream, true);
    if (!self.process_timer_)
      self.process_timer_ = Tcl_CreateTimerHandler(kEventPollMilliseconds, poll_process_io, &self);
  }

  static void delete_process_io(struct ged_subprocess* process, bu_process_io_t stream)
  {
    auto& self = *static_cast<Impl*>(process->gedp->ged_io_data);
    self.process_listeners_.erase({process, stream});
    set_process_stream(process, stream, false);
  }

  static void poll_process_io(ClientData data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.process_timer_ = nullptr;
    const auto listeners = self.process_listeners_;
    for (const auto& entry : listeners) {
      if (!self.process_listeners_.count(entry.first)) continue;
      const auto& listener = entry.second;
      if (listener.stream == BU_PROCESS_STDIN) continue;
      const int descriptor = bu_process_fileno(listener.process->p, listener.stream);
      if (descriptor >= 0 && (bu_process_pending(descriptor) ||
                             bu_process_poll(listener.process->p, nullptr))) {
        listener.callback(listener.data, static_cast<int>(listener.stream));
      }
    }
    if (!self.process_listeners_.empty())
      self.process_timer_ = Tcl_CreateTimerHandler(kEventPollMilliseconds, poll_process_io, &self);
  }

  static int refresh_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                              Tcl_Obj* const values[])
  {
    if (count != 1) {
      Tcl_WrongNumArgs(interpreter, 1, values, "");
      return TCL_ERROR;
    }
    auto& self = *static_cast<Impl*>(data);
    if (self.changed_handler_) self.changed_handler_();
    return TCL_OK;
  }

  static int output_hook_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                                   Tcl_Obj* const values[])
  {
    auto& self = *static_cast<Impl*>(data);
    if (count > 2) {
      Tcl_WrongNumArgs(interpreter, 1, values, "?command?");
      return TCL_ERROR;
    }
    const std::string command = count == 2 ? Tcl_GetString(values[1]) : "";
    Tcl_CmdInfo info;
    if (!command.empty() && !Tcl_GetCommandInfo(interpreter, command.c_str(), &info)) {
      Tcl_SetObjResult(interpreter, Tcl_NewStringObj("output_hook requires an existing Tcl command", -1));
      return TCL_ERROR;
    }
    self.drain_output_hook();
    bu_log_delete_hook(capture_output_hook, &self);
    self.output_hook_command_ = command;
    if (!command.empty()) bu_log_add_hook(capture_output_hook, &self);
    return TCL_OK;
  }

  static int capture_output_hook(void* data, void* message)
  {
    auto& self = *static_cast<Impl*>(data);
    std::lock_guard<std::mutex> lock(self.output_hook_mutex_);
    const std::string text = static_cast<const char*>(message);
    const std::size_t accepted = std::min(text.size(), kMaxResultBytes - self.output_hook_pending_.size());
    self.output_hook_pending_.append(text, 0, accepted);
    if (accepted < text.size()) self.output_hook_truncated_ = true;
    return 0;
  }

  void drain_output_hook()
  {
    if (output_hook_command_.empty() || draining_output_hook_) return;
    std::string output;
    {
      std::lock_guard<std::mutex> lock(output_hook_mutex_);
      output.swap(output_hook_pending_);
      if (output_hook_truncated_) {
        output += "\n[output_hook output truncated]\n";
        output_hook_truncated_ = false;
      }
    }
    if (output.empty()) return;
    draining_output_hook_ = true;
    const Tcl_InterpState state = Tcl_SaveInterpState(interpreter_, TCL_OK);
    const std::string ged_result = ged_ ? bu_vls_cstr(ged_->ged_result_str) : "";
    Tcl_Obj* arguments[] = {Tcl_NewStringObj(output_hook_command_.c_str(), -1),
      Tcl_NewStringObj(output.data(), static_cast<int>(output.size()))};
    for (auto* argument : arguments) Tcl_IncrRefCount(argument);
    if (Tcl_EvalObjv(interpreter_, 2, arguments, TCL_EVAL_GLOBAL) != TCL_OK)
      Tcl_BackgroundError(interpreter_);
    for (auto* argument : arguments) Tcl_DecrRefCount(argument);
    if (ged_) bu_vls_sprintf(ged_->ged_result_str, "%s", ged_result.c_str());
    Tcl_RestoreInterpState(interpreter_, state);
    draining_output_hook_ = false;
  }

  static int search_dispatch(int count, const char** arguments, void*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const Tcl_InterpState state = Tcl_SaveInterpState(self.interpreter_, TCL_OK);
    const std::string output = bu_vls_cstr(self.ged_->ged_result_str);
    std::vector<Tcl_Obj*> values;
    for (int index = 0; index < count; ++index) {
      auto* value = Tcl_NewStringObj(arguments[index], -1);
      Tcl_IncrRefCount(value);
      values.push_back(value);
    }
    const int status = Tcl_EvalObjv(self.interpreter_, count, values.data(), TCL_EVAL_GLOBAL);
    const std::string error = status == TCL_OK ? "" : Tcl_GetStringResult(self.interpreter_);
    for (auto* value : values) Tcl_DecrRefCount(value);
    bu_vls_sprintf(self.ged_->ged_result_str, "%s", output.c_str());
    if (!error.empty()) bu_vls_printf(self.ged_->ged_result_str, "%s\n", error.c_str());
    Tcl_RestoreInterpState(self.interpreter_, state);
    return status == TCL_OK;
  }

  class SearchCallback
  {
  public:
    SearchCallback(Impl& owner, bool enabled) : owner_(owner), enabled_(enabled)
    {
      if (!enabled_) return;
      ged_clbk_get(&callback_, &data_, owner_.ged_, "search", BU_CLBK_DURING);
      ged_clbk_set(owner_.ged_, "search", BU_CLBK_DURING, search_dispatch, &owner_);
    }
    ~SearchCallback()
    {
      if (enabled_) ged_clbk_set(owner_.ged_, "search", BU_CLBK_DURING, callback_, data_);
    }
  private:
    Impl& owner_;
    bool enabled_;
    bu_clbk_t callback_ = nullptr;
    void* data_ = nullptr;
  };

  static int commands_dispatch(ClientData data, Tcl_Interp* interpreter,
                               int count, Tcl_Obj* const arguments[])
  {
    if (count != 1) {
      Tcl_WrongNumArgs(interpreter, 1, arguments, nullptr);
      return TCL_ERROR;
    }
    auto& self = *static_cast<Impl*>(data);
    Tcl_Obj* result = Tcl_NewListObj(0, nullptr);
    for (const auto& command : self.commands_)
      Tcl_ListObjAppendElement(interpreter, result, Tcl_NewStringObj(command->name.c_str(), -1));
    Tcl_SetObjResult(interpreter, result);
    return TCL_OK;
  }

  static char* units_trace(ClientData data, Tcl_Interp* interpreter,
                            const char* name, const char*, int)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->dbip) {
      static char error[] = "No BRL-CAD database is open";
      return error;
    }
    const double value = std::string(name) == "base2local"
      ? self.ged_->dbip->dbi_base2local : self.ged_->dbip->dbi_local2base;
    Tcl_SetVar2Ex(interpreter, name, nullptr, Tcl_NewDoubleObj(value), TCL_GLOBAL_ONLY);
    return nullptr;
  }

  static char* perspective_trace(ClientData data, Tcl_Interp* interpreter,
                                  const char*, const char*, int flags)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ged_ || !self.ged_->ged_gvp) return nullptr;
    if (flags & TCL_TRACE_WRITES) {
      double value = 0.0;
      Tcl_Obj* variable = Tcl_GetVar2Ex(interpreter, "perspective", nullptr, TCL_GLOBAL_ONLY);
      if (!variable || Tcl_GetDoubleFromObj(interpreter, variable, &value) != TCL_OK) {
        static char error[] = "perspective must be a number";
        return error;
      }
      bv_view_set_perspective(self.ged_->ged_gvp, value);
      if (self.changed_handler_) self.changed_handler_();
    } else {
      Tcl_SetVar2Ex(interpreter, "perspective", nullptr,
                    Tcl_NewDoubleObj(self.ged_->ged_gvp->gv_perspective), TCL_GLOBAL_ONLY);
    }
    return nullptr;
  }

  static int puts_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                           Tcl_Obj* const values[])
  {
    auto& self = *static_cast<Impl*>(data);
    int index = 1;
    const bool no_newline = count > index &&
      std::string(Tcl_GetString(values[index])) == "-nonewline";
    if (no_newline) ++index;
    const int remaining = count - index;
    const std::string channel = remaining == 2 ? Tcl_GetString(values[index]) : "stdout";
    if (self.output_handler_ && (remaining == 1 || remaining == 2) &&
        (channel == "stdout" || channel == "stderr")) {
      int length = 0;
      const char* text = Tcl_GetStringFromObj(values[count - 1], &length);
      std::string output(text, static_cast<std::size_t>(length));
      if (!no_newline) output += '\n';
      self.output_handler_(output, channel == "stderr");
      return TCL_OK;
    }
    std::vector<Tcl_Obj*> arguments(values, values + count);
    arguments[0] = Tcl_NewStringObj("::fltkged::tcl_puts", -1);
    Tcl_IncrRefCount(arguments[0]);
    const int status = Tcl_EvalObjv(interpreter, count, arguments.data(), 0);
    Tcl_DecrRefCount(arguments[0]);
    return status;
  }

  static int exit_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                           Tcl_Obj* const values[])
  {
    auto& self = *static_cast<Impl*>(data);
    int exit_status = 0;
    if (count > 2) {
      Tcl_WrongNumArgs(interpreter, 1, values, "?returnCode?");
      return TCL_ERROR;
    }
    if (count == 2 && Tcl_GetIntFromObj(interpreter, values[1], &exit_status) != TCL_OK)
      return TCL_ERROR;
    self.exit_status_ = exit_status;
    return self.execute_ged(interpreter, "quit", 0, nullptr);
  }

  static int ged_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                          Tcl_Obj* const values[])
  {
    auto& self = *static_cast<Impl*>(data);
    if (count < 2) {
      Tcl_WrongNumArgs(interpreter, 1, values, "command ?argument ...?");
      return TCL_ERROR;
    }
    return self.execute_ged(interpreter, Tcl_GetString(values[1]), count - 2,
                            values + 2);
  }

  static int direct_dispatch(ClientData data, Tcl_Interp* interpreter,
                             int count, Tcl_Obj* const values[])
  {
    auto& command = *static_cast<Command*>(data);
    return command.owner->execute_ged(interpreter, command.name, count - 1,
                                      values + 1);
  }

  static int ui_dispatch(ClientData data, Tcl_Interp* interpreter, int count,
                         Tcl_Obj* const values[])
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.ui_handler_) {
      Tcl_SetObjResult(
        interpreter,
        Tcl_NewStringObj("The FLTK user-interface bridge is not connected", -1));
      return TCL_ERROR;
    }
    std::vector<std::string> arguments;
    arguments.reserve(count > 1 ? static_cast<std::size_t>(count - 1) : 0);
    for (int index = 1; index < count; ++index) {
      arguments.emplace_back(Tcl_GetString(values[index]));
    }
    std::string result;
    const int status = self.ui_handler_(arguments, result);
    if (result.size() > kMaxResultBytes) {
      Tcl_SetObjResult(
        interpreter,
        Tcl_NewStringObj("FLTK user-interface result exceeds the configured limit",
                         -1));
      return TCL_ERROR;
    }
    Tcl_SetObjResult(interpreter,
                     Tcl_NewStringObj(result.data(), static_cast<int>(result.size())));
    return status == 0 ? TCL_OK : TCL_ERROR;
  }

  int execute_ged(Tcl_Interp* interpreter, const std::string& command,
                  int argument_count, Tcl_Obj* const values[])
  {
    if (!ged_) {
      Tcl_SetObjResult(interpreter,
                       Tcl_NewStringObj("No BRL-CAD database is open", -1));
      return TCL_ERROR;
    }
    if (argument_count >= kMaxGedArguments) {
      Tcl_SetObjResult(interpreter, Tcl_NewStringObj("GED argument limit exceeded", -1));
      return TCL_ERROR;
    }
    std::vector<const char*> arguments;
    arguments.reserve(static_cast<std::size_t>(argument_count) + 2);
    arguments.push_back(command.c_str());
    for (int index = 0; index < argument_count; ++index) {
      int length = 0;
      const char* value = Tcl_GetStringFromObj(values[index], &length);
      if (std::string(value).size() != static_cast<std::size_t>(length)) {
        Tcl_SetObjResult(interpreter, Tcl_NewStringObj("GED arguments cannot contain null bytes", -1));
        return TCL_ERROR;
      }
      arguments.push_back(value);
    }
    std::vector<std::string> expanded;
    constexpr const char* glob_commands[] = {
      "draw", "e", "erase", "d", "B", "blast", "kill", "killall", "killtree",
      "l", "g", "group", "r", "comb", "c", "hide", "unhide"
    };
    bool expand = false;
    for (const char* name : glob_commands) if (command == name) expand = true;
    if (expand && ged_->dbip) {
      for (std::size_t index = 1; index < arguments.size(); ++index) {
        const std::string argument = arguments[index];
        if (argument.find_first_of("*?[") == std::string::npos) {
          expanded.push_back(argument);
          continue;
        }
        const char* glob_arguments[] = {"expand", argument.c_str(), nullptr};
        bu_vls_trunc(ged_->ged_result_str, 0);
        if (ged_exec_expand(ged_, 2, glob_arguments) != BRLCAD_OK) {
          Tcl_SetObjResult(interpreter, Tcl_NewStringObj(bu_vls_cstr(ged_->ged_result_str), -1));
          return TCL_ERROR;
        }
        int match_count = 0;
        const char** matches = nullptr;
        if (Tcl_SplitList(interpreter, bu_vls_cstr(ged_->ged_result_str),
                          &match_count, &matches) != TCL_OK) return TCL_ERROR;
        if (expanded.size() + static_cast<std::size_t>(match_count) >= kMaxGedArguments) {
          Tcl_Free(reinterpret_cast<char*>(matches));
          Tcl_SetObjResult(interpreter, Tcl_NewStringObj("Expanded GED argument limit exceeded", -1));
          return TCL_ERROR;
        }
        for (int match = 0; match < match_count; ++match) expanded.emplace_back(matches[match]);
        Tcl_Free(reinterpret_cast<char*>(matches));
        if (!match_count) expanded.push_back(argument);
      }
      arguments.resize(1);
      for (const auto& argument : expanded) arguments.push_back(argument.c_str());
    }
    std::string origin;
    std::string scale;
    if (command == "make" && arguments.size() == 3 && ged_->ged_gvp) {
      const bview* view = ged_->ged_gvp;
      std::ostringstream coordinates;
      coordinates << std::setprecision(17) << -view->gv_center[MDX] << ' '
                  << -view->gv_center[MDY] << ' ' << -view->gv_center[MDZ];
      origin = coordinates.str();
      std::ostringstream extent;
      extent << std::setprecision(17) << view->gv_size;
      scale = extent.str();
      arguments.insert(arguments.begin() + 1, {"-o", origin.c_str(), "-s", scale.c_str()});
    }
    if (command == "shaded_mode" && arguments.size() == 3 &&
        (std::string(arguments[1]) == "-a" || std::string(arguments[1]) == "-auto"))
      arguments.erase(arguments.begin() + 1);
    // MGED accepts an explicit affirmative response when creating a database.
    if (command == "opendb" && arguments.size() == 3 &&
        std::string(arguments.back()) == "y") arguments.pop_back();
    arguments.push_back(nullptr);
    bu_vls_trunc(ged_->ged_result_str, 0);
    SearchCallback search_callback(*this, command == "search");
    const int status =
      ged_edit_session_exec(ged_, static_cast<int>(arguments.size() - 1), arguments.data());
    if (command == "echo" && status == BRLCAD_OK) {
      bu_log("%s", bu_vls_cstr(ged_->ged_result_str));
      bu_vls_trunc(ged_->ged_result_str, 0);
    }
    drain_output_hook();
    command_status_ = status;
    const std::size_t size = bu_vls_strlen(ged_->ged_result_str);
    if (size > kMaxResultBytes ||
        size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      bu_vls_trunc(ged_->ged_result_str, 0);
      Tcl_SetObjResult(interpreter,
                       Tcl_NewStringObj("GED result exceeds the configured limit", -1));
      return TCL_ERROR;
    }
    Tcl_SetObjResult(
      interpreter,
      Tcl_NewStringObj(bu_vls_cstr(ged_->ged_result_str), static_cast<int>(size)));
    bu_vls_trunc(ged_->ged_result_str, 0);
    if (status & GED_EXIT) {
      exit_requested_ = true;
      Tcl_CancelEval(interpreter, nullptr, nullptr, TCL_CANCEL_UNWIND);
      // Consume the cancellation before this evaluation returns, so reopening
      // a window cannot deliver the old request to the next command.
      return Tcl_AsyncInvoke(interpreter, TCL_ERROR);
    }
    if ((status & ~(GED_EXIT | GED_HELP)) == BRLCAD_OK) {
      if (changed_handler_) {
        changed_handler_();
      }
      return TCL_OK;
    }
    return TCL_ERROR;
  }

  void copy_result()
  {
    int size = 0;
    const char* value = Tcl_GetStringFromObj(Tcl_GetObjResult(interpreter_), &size);
    if (!value || size <= 0) {
      result_.clear();
      return;
    }
    if (static_cast<std::size_t>(size) > kMaxResultBytes) {
      result_.assign(value, kMaxResultBytes);
      result_.append("\n[Tcl result truncated by application limit]");
      last_error_ = "Tcl result exceeds the configured size limit";
      return;
    }
    result_.assign(value, static_cast<std::size_t>(size));
  }

  bool set_limits()
  {
    const Tcl_InterpState state = Tcl_SaveInterpState(interpreter_, TCL_OK);
    int commands = 0;
    const bool valid = Tcl_EvalEx(interpreter_, "::tcl::info::cmdcount", -1, TCL_EVAL_GLOBAL) == TCL_OK &&
      Tcl_GetIntFromObj(interpreter_, Tcl_GetObjResult(interpreter_), &commands) == TCL_OK;
    Tcl_RestoreInterpState(interpreter_, state);
    if (!valid) {
      last_error_ = "Unable to read Tcl's command count for the evaluation limit";
      return false;
    }
    Tcl_Time deadline;
    Tcl_GetTime(&deadline);
    constexpr long kEvaluationSeconds = 30;
    deadline.sec += kEvaluationSeconds;
    Tcl_LimitSetTime(interpreter_, &deadline);
    Tcl_LimitTypeSet(interpreter_, TCL_LIMIT_TIME);
    constexpr int kEvaluationCommands = 1000000;
    // Tcl's limit is an absolute interpreter count, rather than a per-script budget.
    const int limit = commands > std::numeric_limits<int>::max() - kEvaluationCommands
      ? std::numeric_limits<int>::max() : commands + kEvaluationCommands;
    Tcl_LimitSetCommands(interpreter_, limit);
    Tcl_LimitTypeSet(interpreter_, TCL_LIMIT_COMMANDS);
    return true;
  }

  void clear_limits()
  {
    Tcl_LimitTypeReset(interpreter_, TCL_LIMIT_TIME);
    Tcl_LimitTypeReset(interpreter_, TCL_LIMIT_COMMANDS);
  }

  void start_tk_events_if_loaded()
  {
    if (tk_active_ || !interpreter_ ||
        !Tcl_PkgPresent(interpreter_, "Tk", nullptr, 0)) {
      return;
    }
    tk_active_ = true;
    event_processing(true);
  }

  static void service_tcl_events(void* data)
  {
    auto* self = static_cast<Impl*>(data);
    if (!self || !self->events_enabled_ || !self->interpreter_) {
      return;
    }
    if (!self->evaluating_ && !self->servicing_events_) {
      self->servicing_events_ = true;
      if (!self->set_limits()) {
        self->servicing_events_ = false;
        self->event_processing(false);
        if (self->output_handler_) self->output_handler_(self->last_error_ + "\n", true);
        return;
      }
      for (int count = 0; count < kMaxTclEventsPerTick; ++count) {
        if (!Tcl_DoOneEvent(TCL_ALL_EVENTS | TCL_DONT_WAIT)) break;
      }
      self->clear_limits();
      self->drain_output_hook();
      self->servicing_events_ = false;
    }
    Fl::repeat_timeout(kTclEventIntervalSeconds, service_tcl_events, self);
  }

  struct FltkEvent { Tcl_Event event; Impl* owner; };

  static void setup_fltk_events(ClientData data, int flags)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!(flags & TCL_WINDOW_EVENTS) || self.fltk_check_active_) return;
    Tcl_Time delay{0, Fl::ready() ? 0 : kEventPollMilliseconds * 1000};
    Tcl_SetMaxBlockTime(&delay);
  }

  static void check_fltk_events(ClientData data, int flags)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!(flags & TCL_WINDOW_EVENTS) || self.fltk_check_active_ ||
        self.fltk_event_queued_) return;
    if (!Fl::ready() && std::chrono::steady_clock::now() < self.next_fltk_poll_)
      return;
    auto* event = reinterpret_cast<FltkEvent*>(Tcl_Alloc(sizeof(FltkEvent)));
    event->event.proc = service_fltk_event;
    event->owner = &self;
    self.fltk_event_queued_ = true;
    Tcl_QueueEvent(&event->event, TCL_QUEUE_TAIL);
  }

  static int service_fltk_event(Tcl_Event* event, int flags)
  {
    if (!(flags & TCL_WINDOW_EVENTS)) return 0;
    auto& self = *reinterpret_cast<FltkEvent*>(event)->owner;
    self.fltk_event_queued_ = false;
    self.next_fltk_poll_ = std::chrono::steady_clock::now() +
      std::chrono::milliseconds(kEventPollMilliseconds);
    self.fltk_check_active_ = true;
    Fl::check();
    self.fltk_check_active_ = false;
    return 1;
  }

  static int delete_fltk_event(Tcl_Event* event, ClientData data)
  {
    return event->proc == service_fltk_event &&
      reinterpret_cast<FltkEvent*>(event)->owner == data;
  }

  Tcl_Interp* interpreter_ = nullptr;
  std::string output_hook_command_;
  std::mutex output_hook_mutex_;
  std::string output_hook_pending_;
  bool draining_output_hook_ = false;
  bool output_hook_truncated_ = false;
  bool owns_process_io_ = false;
  Tcl_TimerToken process_timer_ = nullptr;
  std::map<ProcessKey, ProcessListener> process_listeners_;
  bool fltk_event_queued_ = false;
  bool fltk_check_active_ = false;
  std::chrono::steady_clock::time_point next_fltk_poll_;
  std::vector<std::unique_ptr<Command>> commands_;
#endif
};

Fl_Ged_Tcl::Fl_Ged_Tcl() : impl_(std::make_unique<Impl>()) {}
Fl_Ged_Tcl::~Fl_Ged_Tcl() = default;
bool Fl_Ged_Tcl::available() const noexcept { return impl_->available_; }
void Fl_Ged_Tcl::ged(struct ged* value) { impl_->bind_ged(value); }
int Fl_Ged_Tcl::evaluate(const std::string& script)
{
  return impl_->evaluate(script);
}
int Fl_Ged_Tcl::source(const std::string& path) { return impl_->source(path); }
void Fl_Ged_Tcl::interrupt() { impl_->interrupt(); }
bool Fl_Ged_Tcl::command_complete(const std::string& script) const
{
  return impl_->command_complete(script);
}
int Fl_Ged_Tcl::command_status() const noexcept { return impl_->command_status_; }
bool Fl_Ged_Tcl::exit_requested() const noexcept { return impl_->exit_requested_; }
int Fl_Ged_Tcl::exit_status() const noexcept { return impl_->exit_status_; }
void Fl_Ged_Tcl::output_handler(OutputHandler handler)
{
  impl_->output_handler_ = std::move(handler);
}
void Fl_Ged_Tcl::event_processing(bool enabled) { impl_->event_processing(enabled); }
bool Fl_Ged_Tcl::enable_legacy_tk() { return impl_->enable_legacy_tk(); }
bool Fl_Ged_Tcl::legacy_tk_active() const noexcept { return impl_->tk_active_; }
const std::string& Fl_Ged_Tcl::result() const noexcept { return impl_->result_; }
const std::string& Fl_Ged_Tcl::last_error() const noexcept
{
  return impl_->last_error_;
}
void Fl_Ged_Tcl::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler_ = std::move(handler);
}
void Fl_Ged_Tcl::ui_handler(UiHandler handler)
{
  impl_->ui_handler_ = std::move(handler);
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
