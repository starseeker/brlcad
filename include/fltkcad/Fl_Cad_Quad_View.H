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

#ifndef FLTKCAD_FL_CAD_QUAD_VIEW_H
#define FLTKCAD_FL_CAD_QUAD_VIEW_H

#include "fltkcad/defines.h"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Group.H>

#include <array>
#include <cstddef>
#include <functional>

class FLTKCAD_EXPORT Fl_Cad_Quad_View : public Fl_Group
{
public:
  enum class Layout
  {
    Single,
    Quad
  };

  using ActiveHandler = std::function<void(std::size_t index)>;
  using ChangedHandler = std::function<void(std::size_t index)>;
  using InputHandler =
    std::function<bool(std::size_t index, const Fl_Cad_View::InputEvent& event)>;
  using ErrorHandler =
    std::function<void(std::size_t index, const std::string& message)>;

  Fl_Cad_Quad_View(int x, int y, int width, int height,
                   Fl_Cad_View::Backend backend =
                     Fl_Cad_View::Backend::Auto,
                   const char* label = nullptr);
  ~Fl_Cad_Quad_View() override;

  Fl_Cad_Quad_View(const Fl_Cad_Quad_View&) = delete;
  Fl_Cad_Quad_View& operator=(const Fl_Cad_Quad_View&) = delete;

  Fl_Cad_View& active_view() noexcept;
  const Fl_Cad_View& active_view() const noexcept;
  Fl_Cad_View& view(std::size_t index);
  const Fl_Cad_View& view(std::size_t index) const;
  std::array<Fl_Cad_View*, 4> views() noexcept;

  std::size_t active_index() const noexcept;
  void active_index(std::size_t index);

  Layout layout() const noexcept;
  void layout(Layout value);

  void active_handler(ActiveHandler handler);
  void input_handler(InputHandler handler);
  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);

  void invalidate();
  void resize(int x, int y, int width, int height) override;

private:
  class View;
  void activate(std::size_t index, bool changed);
  void arrange();

  std::array<View*, 4> views_{{nullptr, nullptr, nullptr, nullptr}};
  std::size_t active_index_ = 0;
  Layout layout_ = Layout::Single;
  ActiveHandler active_handler_;
  InputHandler input_handler_;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
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
