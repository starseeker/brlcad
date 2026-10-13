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
#include "bu/app.h"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

#include <cstdlib>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

extern "C" {
#include "dm.h"
#include "dm/faceplate.h"
}

namespace {

constexpr double kCoordinateTolerance = 1.0;

} // namespace

int
main(int argc, char** argv)
{
  bu_setprogname(argv[0]);
  if (std::getenv("DISPLAY") == nullptr || argc != 2) {
    return 77;
  }
  const auto backend = std::string(argv[1]) == "opengl"
                         ? Fl_Cad_View::Backend::SystemOpenGL
                         : Fl_Cad_View::Backend::OSMesa;

  Fl_Double_Window window(320, 240);
  Fl_Cad_View view(0, 0, window.w(), window.h(), backend);
  bool text_bounds_checked = backend != Fl_Cad_View::Backend::SystemOpenGL;
  view.draw_handler([&text_bounds_checked](struct bview* bv) {
    auto* dmp = static_cast<struct dm*>(bv->dmp);
    dm_hud_begin(dmp);
    dm_set_fg(dmp, 255, 255, 0, 1, 1.0);
    dm_draw_line_2d(dmp, -0.9, 0.0, 0.9, 0.0);
    dm_draw_line_2d(dmp, 0.0, -0.9, 0.0, 0.9);
    if (!text_bounds_checked) {
      vect2d_t minimum = V2INIT_ZERO;
      vect2d_t maximum = V2INIT_ZERO;
      const double center_x = dm_get_width(dmp) * 0.5;
      const double center_y = dm_get_height(dmp) * 0.5;
      if (dm_string_bbox_2d(dmp, &minimum, &maximum, "label", 0.0, 0.0,
                            12, 0) != BRLCAD_OK ||
          std::abs(minimum[0] - center_x) > kCoordinateTolerance ||
          std::abs(maximum[1] - center_y) > kCoordinateTolerance ||
          maximum[0] <= minimum[0] || maximum[1] <= minimum[1]) {
        throw std::runtime_error("text bounds are not screen coordinates");
      }
      text_bounds_checked = true;
    }
    dm_hud_end(dmp);
  });
  window.end();
  window.show();
  view.invalidate();
  const auto deadline = std::chrono::steady_clock::now() +
                        std::chrono::seconds(2);
  while (!view.initialized() && std::chrono::steady_clock::now() < deadline) {
    Fl::wait(0.01);
  }

  if (!view.initialized() && !view.initialize()) {
    std::cerr << view.last_error() << '\n';
    return EXIT_FAILURE;
  }
  auto image = view.capture();
  if (!image || image->width <= 0 || image->height <= 0 ||
      image->rgba.size() !=
        static_cast<std::size_t>(image->width) * image->height * 4U) {
    std::cerr << view.last_error() << '\n';
    return EXIT_FAILURE;
  }
  bool nonzero_rgb = false;
  for (std::size_t offset = 0; offset + 3 < image->rgba.size(); offset += 4) {
    if (image->rgba[offset] != 0 || image->rgba[offset + 1] != 0 ||
        image->rgba[offset + 2] != 0) {
      nonzero_rgb = true;
      break;
    }
  }
  if (!nonzero_rgb) {
    std::cerr << "captured display is blank\n";
    return EXIT_FAILURE;
  }
  if (!text_bounds_checked) {
    std::cerr << "text bounds were not checked\n";
    return EXIT_FAILURE;
  }

  struct dm_faceplate_controls menu;
  dm_faceplate_controls_init(&menu, DM_FACEPLATE_MENU);
  menu.item = [](void*, std::size_t, struct dm_faceplate_item*) -> int {
    throw std::runtime_error("intentional faceplate provider failure");
  };
  struct dm_faceplate_gui gui{&menu, nullptr};
  view.faceplate_gui(&gui);
  if (view.capture() || view.last_error().find("Faceplate drawing callback failed") == std::string::npos) {
    std::cerr << "faceplate provider exception was not contained\n";
    return EXIT_FAILURE;
  }
  view.faceplate_gui(nullptr);
  if (!view.capture()) {
    std::cerr << "display did not recover after unregistering failed faceplate\n";
    return EXIT_FAILURE;
  }

  view.error_handler([](const std::string&) {
    throw std::runtime_error("intentional error-handler failure");
  });
  view.draw_handler([](struct bview*) {
    throw std::runtime_error("intentional drawing failure");
  });
  view.invalidate();
  Fl::wait(0.05);
  if (view.last_error().find("Display drawing callback failed") == std::string::npos) {
    std::cerr << "drawing callback exception was not reported\n";
    return EXIT_FAILURE;
  }

  auto limits = view.limits();
  limits.max_framebuffer_width = 1;
  view.limits(limits);
  if (view.capture() || view.last_error().find("dimensions") == std::string::npos) {
    std::cerr << "framebuffer limit failure was not reported\n";
    return EXIT_FAILURE;
  }
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
