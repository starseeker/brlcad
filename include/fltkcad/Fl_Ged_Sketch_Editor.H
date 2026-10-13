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

#ifndef FLTKCAD_FL_GED_SKETCH_EDITOR_H
#define FLTKCAD_FL_GED_SKETCH_EDITOR_H

#include "fltkcad/defines.h"

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl_Double_Window.H>

#include <functional>
#include <memory>
#include <string>
#include <vector>

struct ged;

class FLTKCAD_EXPORT Fl_Ged_Sketch_Editor : public Fl_Double_Window
{
public:
  enum class VertexConstraint
  {
    Coincident,
    Horizontal,
    Vertical,
    Distance
  };

  using ChangedHandler = std::function<void()>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  explicit Fl_Ged_Sketch_Editor(
    Fl_Cad_View::Backend backend = Fl_Cad_View::Backend::Auto);
  ~Fl_Ged_Sketch_Editor() override;

  Fl_Ged_Sketch_Editor(const Fl_Ged_Sketch_Editor&) = delete;
  Fl_Ged_Sketch_Editor& operator=(const Fl_Ged_Sketch_Editor&) = delete;

  bool open(struct ged* gedp, const std::string& object_name);
  bool save();
  bool undo();
  bool move_vertex(std::size_t index, double u, double v);
  bool move_vertices(const std::vector<std::size_t>& indices, double delta_u,
                     double delta_v);
  bool constrain_vertices(VertexConstraint constraint, std::size_t fixed,
                          std::size_t moving, double distance = 0.0);
  bool split_segment(std::size_t index, double parameter);
  bool reverse_segment(std::size_t index);
  bool set_arc_radius(std::size_t index, double radius);
  bool toggle_arc_side(std::size_t index);
  bool set_tangency(std::size_t arc, std::size_t adjacent,
                    double angle_radians = 0.0);
  bool set_nurb_knots(std::size_t index,
                      const std::vector<double>& knots);
  bool set_nurb_weights(std::size_t index,
                        const std::vector<double>& weights);

  bool dirty() const noexcept;
  const std::string& object_name() const noexcept;
  const std::string& last_error() const noexcept;
  void changed_handler(ChangedHandler handler);
  void error_handler(ErrorHandler handler);
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
