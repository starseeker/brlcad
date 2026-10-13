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

#ifndef FLTKCAD_COMMAND_HISTORY_H
#define FLTKCAD_COMMAND_HISTORY_H

#include "fltkcad/defines.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace fltkcad {

class FLTKCAD_EXPORT CommandHistory
{
public:
  explicit CommandHistory(std::size_t limit = 1000);

  void add(std::string command);
  void clear();
  void limit(std::size_t value);
  std::size_t limit() const noexcept;

  std::optional<std::string> previous(const std::string& current_input);
  std::optional<std::string> next();
  void reset_navigation();

  bool consolidate(std::size_t begin, std::size_t end);
  bool consolidate(std::size_t begin, std::size_t end, std::string command);
  std::size_t size() const noexcept;
  const std::string& at(std::size_t index) const;

private:
  std::size_t limit_;
  std::vector<std::string> entries_;
  std::size_t position_ = 0;
  std::string draft_;
};

} // namespace fltkcad

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
