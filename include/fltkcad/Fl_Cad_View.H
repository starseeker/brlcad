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

#ifndef FLTKCAD_FL_CAD_VIEW_H
#define FLTKCAD_FL_CAD_VIEW_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

struct bview;
struct dm;
struct fb;
struct dm_faceplate_gui;

struct ged;

class FLTKCAD_EXPORT Fl_Cad_View : public Fl_Group
{
public:
  enum class Backend
  {
    Auto,
    SystemOpenGL,
    OSMesa
  };

  struct Image
  {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgba;
  };

  struct Color
  {
    unsigned char red = 0;
    unsigned char green = 0;
    unsigned char blue = 0;
  };

  struct Limits
  {
    int max_framebuffer_width = 4096;
    int max_framebuffer_height = 4096;
    std::size_t max_framebuffer_pixels = 16U * 1024U * 1024U;
  };

  struct InputEvent
  {
    int type = 0;
    int x = 0;
    int y = 0;
    int button = 0;
    int buttons = 0;
    int state = 0;
    int wheel_x = 0;
    int wheel_y = 0;
    int key = 0;
  };

  using DrawHandler = std::function<void(struct bview* view)>;
  using InputHandler = std::function<bool(const InputEvent& event)>;
  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Cad_View(int x, int y, int width, int height,
              Backend preference = Backend::Auto,
              const char* label = nullptr);
  ~Fl_Cad_View() override;

  Fl_Cad_View(const Fl_Cad_View&) = delete;
  Fl_Cad_View& operator=(const Fl_Cad_View&) = delete;

  bool initialize();
  bool initialized() const noexcept;
  Backend backend() const noexcept;
  Backend preference() const noexcept;

  /** Borrow a GED context to display its active working edit. */
  void ged(struct ged* context);
  struct bview* view() noexcept;
  const struct bview* view() const noexcept;
  void view(struct bview* borrowed_view);
  struct dm* display_manager() noexcept;
  struct fb* framebuffer() noexcept;

  // Borrowed description; unregister before destroying its callback data.
  void faceplate_gui(const struct dm_faceplate_gui* gui);
  const struct dm_faceplate_gui* faceplate_gui() const noexcept;

  void draw_handler(DrawHandler handler);
  void input_handler(InputHandler handler);
  bool process_input(const InputEvent& event);
  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);

  void default_bindings(bool enabled) noexcept;
  bool default_bindings() const noexcept;

  void limits(Limits value);
  Limits limits() const noexcept;
  void background(Color first, Color second);
  Color background_first() const noexcept;
  Color background_second() const noexcept;
  const std::string& last_error() const noexcept;

  void invalidate();
  std::optional<Image> capture();

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
