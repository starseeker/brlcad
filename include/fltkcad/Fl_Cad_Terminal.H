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

#ifndef FLTKCAD_FL_CAD_TERMINAL_H
#define FLTKCAD_FL_CAD_TERMINAL_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

class Fl_Terminal;

class FLTKCAD_EXPORT Fl_Cad_Terminal : public Fl_Group
{
public:
  struct Completion
  {
    std::size_t begin = 0;
    std::size_t end = 0;
    std::vector<std::string> candidates;
  };

  struct Limits
  {
    std::size_t command_history = 1000;
    int scrollback_lines = 10000;
    std::size_t max_input_bytes = 64U * 1024U;
    std::size_t max_append_bytes = 4U * 1024U * 1024U;
    std::size_t max_prompt_bytes = 4U * 1024U;
    std::size_t max_completion_candidates = 200;
  };

  using CommandHandler = std::function<void(const std::string& command)>;
  using CompletionHandler =
    std::function<Completion(const std::string& input, std::size_t cursor)>;
  using InterruptHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Cad_Terminal(int x, int y, int width, int height,
                  const char* label = nullptr);
  ~Fl_Cad_Terminal() override;

  Fl_Cad_Terminal(const Fl_Cad_Terminal&) = delete;
  Fl_Cad_Terminal& operator=(const Fl_Cad_Terminal&) = delete;

  void append(const std::string& text);
  void clear();
  /// Current viewport width in character cells, independent of scrollback width.
  int columns() const noexcept;
  /// Change the default text and background colors, including existing text.
  void colors(Fl_Color foreground, Fl_Color background);

  void prompt(std::string text);
  const std::string& prompt() const noexcept;

  void input(std::string text);
  std::string input() const;
  void submit();
  void complete();

  void command_handler(CommandHandler handler);
  void completion_handler(CompletionHandler handler);
  void interrupt_handler(InterruptHandler handler);
  void error_handler(ErrorHandler handler);

  void limits(Limits value);
  Limits limits() const noexcept;
  const std::string& last_error() const noexcept;

  std::size_t history_size() const;
  std::string history_at(std::size_t index) const;
  bool consolidate_history(std::size_t begin, std::size_t end);
  bool consolidate_history(std::size_t begin, std::size_t end, const std::string& command);

  Fl_Terminal& output() noexcept;
  const Fl_Terminal& output() const noexcept;

  void resize(int x, int y, int width, int height) override;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
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
