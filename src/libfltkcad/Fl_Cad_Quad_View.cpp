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

#include "fltkcad/Fl_Cad_Quad_View.H"

#include <FL/Fl.H>

#include <algorithm>
#include <stdexcept>
#include <utility>

extern "C" {
#include "bv.h"
}

namespace {

constexpr std::size_t kViewCount = 4;
constexpr int kDivider = 2;

const std::array<std::array<double, 3>, kViewCount> kDefaultAet{{
  {{270.0, 90.0, 0.0}},
  {{35.0, 25.0, 0.0}},
  {{0.0, 0.0, 0.0}},
  {{90.0, 0.0, 0.0}},
}};

} // namespace

class Fl_Cad_Quad_View::View : public Fl_Cad_View
{
public:
  View(Fl_Cad_Quad_View& owner, std::size_t index, int x, int y, int width,
       int height, Backend backend)
    : Fl_Cad_View(x, y, width, height, backend), owner_(owner), index_(index)
  {
  }

  int handle(int event) override
  {
    if (event == FL_PUSH || event == FL_FOCUS) {
      owner_.activate(index_, false);
    }
    return Fl_Cad_View::handle(event);
  }

private:
  Fl_Cad_Quad_View& owner_;
  std::size_t index_;
};

Fl_Cad_Quad_View::Fl_Cad_Quad_View(int x, int y, int width, int height,
                                   Fl_Cad_View::Backend backend,
                                   const char* label)
  : Fl_Group(x, y, width, height, label)
{
  begin();
  for (std::size_t index = 0; index < views_.size(); ++index) {
    views_[index] = new View(*this, index, x, y, width, height, backend);
    auto* bview = views_[index]->view();
    if (bview) {
      VSET(bview->gv_aet, kDefaultAet[index][0], kDefaultAet[index][1],
           kDefaultAet[index][2]);
      bv_mat_aet(bview);
      bv_update(bview);
    }
    views_[index]->changed_handler([this, index]() {
      activate(index, true);
    });
    views_[index]->input_handler(
      [this, index](const Fl_Cad_View::InputEvent& event) {
        return input_handler_ ? input_handler_(index, event) : false;
      });
    views_[index]->error_handler(
      [this, index](const std::string& message) {
        if (error_handler_) {
          error_handler_(index, message);
        }
      });
  }
  end();
  resizable(views_[0]);
  arrange();
}

Fl_Cad_Quad_View::~Fl_Cad_Quad_View() = default;

Fl_Cad_View& Fl_Cad_Quad_View::active_view() noexcept
{
  return *views_[active_index_];
}

const Fl_Cad_View& Fl_Cad_Quad_View::active_view() const noexcept
{
  return *views_[active_index_];
}

Fl_Cad_View& Fl_Cad_Quad_View::view(std::size_t index)
{
  if (index >= views_.size()) {
    throw std::out_of_range("CAD view index is out of range");
  }
  return *views_[index];
}

const Fl_Cad_View& Fl_Cad_Quad_View::view(std::size_t index) const
{
  if (index >= views_.size()) {
    throw std::out_of_range("CAD view index is out of range");
  }
  return *views_[index];
}

std::array<Fl_Cad_View*, 4> Fl_Cad_Quad_View::views() noexcept
{
  return {{views_[0], views_[1], views_[2], views_[3]}};
}

std::size_t Fl_Cad_Quad_View::active_index() const noexcept
{
  return active_index_;
}

void Fl_Cad_Quad_View::active_index(std::size_t index)
{
  if (index >= views_.size()) {
    return;
  }
  activate(index, false);
}

Fl_Cad_Quad_View::Layout Fl_Cad_Quad_View::layout() const noexcept
{
  return layout_;
}

void Fl_Cad_Quad_View::layout(Layout value)
{
  if (layout_ == value) {
    return;
  }
  layout_ = value;
  arrange();
  redraw();
}

void Fl_Cad_Quad_View::active_handler(ActiveHandler handler)
{
  active_handler_ = std::move(handler);
}

void Fl_Cad_Quad_View::input_handler(InputHandler handler)
{
  input_handler_ = std::move(handler);
}

void Fl_Cad_Quad_View::changed_handler(ChangedHandler handler)
{
  changed_handler_ = std::move(handler);
}

void Fl_Cad_Quad_View::error_handler(ErrorHandler handler)
{
  error_handler_ = std::move(handler);
}

void Fl_Cad_Quad_View::invalidate()
{
  for (auto& view : views_) {
    if (view->visible()) {
      view->invalidate();
    }
  }
}

void Fl_Cad_Quad_View::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  arrange();
}

void Fl_Cad_Quad_View::activate(std::size_t index, bool changed)
{
  if (index >= views_.size()) {
    return;
  }
  const bool active_changed = active_index_ != index;
  active_index_ = index;
  if (active_changed && layout_ == Layout::Single) {
    arrange();
    redraw();
  }
  if (active_changed && active_handler_) {
    active_handler_(index);
  }
  if (changed && changed_handler_) {
    changed_handler_(index);
  }
}

void Fl_Cad_Quad_View::arrange()
{
  if (layout_ == Layout::Single) {
    for (std::size_t index = 0; index < views_.size(); ++index) {
      if (index == active_index_) {
        views_[index]->resize(x(), y(), w(), h());
        views_[index]->show();
      } else {
        views_[index]->hide();
      }
    }
    return;
  }

  const int left_width = std::max(1, (w() - kDivider) / 2);
  const int right_width = std::max(1, w() - kDivider - left_width);
  const int top_height = std::max(1, (h() - kDivider) / 2);
  const int bottom_height = std::max(1, h() - kDivider - top_height);
  views_[0]->resize(x(), y(), left_width, top_height);
  views_[1]->resize(x() + left_width + kDivider, y(), right_width,
                    top_height);
  views_[2]->resize(x(), y() + top_height + kDivider, left_width,
                    bottom_height);
  views_[3]->resize(x() + left_width + kDivider,
                    y() + top_height + kDivider, right_width, bottom_height);
  for (auto& view : views_) {
    view->show();
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
