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

#include "fltkcad/CommandHistory.H"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace fltkcad {

CommandHistory::CommandHistory(std::size_t limit) : limit_(limit)
{
}

void
CommandHistory::add(std::string command)
{
  if (command.empty()) {
    reset_navigation();
    return;
  }

  if (entries_.empty() || entries_.back() != command) {
    entries_.push_back(std::move(command));
    if (entries_.size() > limit_) {
      const auto excess = entries_.size() - limit_;
      entries_.erase(entries_.begin(), entries_.begin() + excess);
    }
  }
  reset_navigation();
}

void
CommandHistory::clear()
{
  entries_.clear();
  reset_navigation();
}

void
CommandHistory::limit(std::size_t value)
{
  limit_ = value;
  if (entries_.size() > limit_) {
    const auto excess = entries_.size() - limit_;
    entries_.erase(entries_.begin(), entries_.begin() + excess);
  }
  reset_navigation();
}

std::size_t
CommandHistory::limit() const noexcept
{
  return limit_;
}

std::optional<std::string>
CommandHistory::previous(const std::string& current_input)
{
  if (entries_.empty()) {
    return std::nullopt;
  }
  if (position_ == entries_.size()) {
    draft_ = current_input;
  }
  if (position_ > 0) {
    --position_;
  }
  return entries_[position_];
}

std::optional<std::string>
CommandHistory::next()
{
  if (position_ >= entries_.size()) {
    return std::nullopt;
  }
  ++position_;
  return position_ == entries_.size() ? draft_ : entries_[position_];
}

void
CommandHistory::reset_navigation()
{
  position_ = entries_.size();
  draft_.clear();
}

bool
CommandHistory::consolidate(std::size_t begin, std::size_t end)
{
  if (begin >= end || end > entries_.size()) {
    return false;
  }

  std::string combined;
  for (std::size_t index = begin; index < end; ++index) {
    if (!combined.empty()) {
      combined.push_back(' ');
    }
    combined += entries_[index];
  }

  return consolidate(begin, end, std::move(combined));
}

bool
CommandHistory::consolidate(std::size_t begin, std::size_t end, std::string command)
{
  if (begin >= end || end > entries_.size()) return false;
  entries_.erase(entries_.begin() + begin, entries_.begin() + end);
  entries_.insert(entries_.begin() + begin, std::move(command));
  reset_navigation();
  return true;
}

std::size_t
CommandHistory::size() const noexcept
{
  return entries_.size();
}

const std::string&
CommandHistory::at(std::size_t index) const
{
  if (index >= entries_.size()) {
    throw std::out_of_range("terminal history index is out of range");
  }
  return entries_[index];
}

} // namespace fltkcad

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
