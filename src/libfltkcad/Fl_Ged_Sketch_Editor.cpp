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
#include "vmath.h"

#include "fltkcad/Fl_Ged_Sketch_Editor.H"

#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Input.H>
#include <FL/fl_ask.H>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "bg.h"
#include "bu/vls.h"
#include "bv/util.h"
#include "bv/vlist.h"
#include "dm.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/edit.h"
#include "rt/functab.h"
#include "rt/primitives/sketch.h"
#include "rt/rt_ecmds.h"
}

namespace {

constexpr int kInitialWidth = 1100;
constexpr int kInitialHeight = 700;
constexpr int kPanelWidth = 340;
constexpr int kMargin = 8;
constexpr int kRowHeight = 27;
constexpr int kButtonHeight = 29;
constexpr std::size_t kMaxVertices = 250000;
constexpr std::size_t kMaxSegments = 250000;
constexpr std::size_t kMaxLogBytes = 64U * 1024U;
constexpr int kPreviewRed = 255;
constexpr int kPreviewGreen = 215;
constexpr int kPreviewBlue = 0;
constexpr double kPickRadiusPixels = 12.0;

bool parse_number(const char* text, double& value)
{
  if (!text || !*text) {
    return false;
  }
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text, &end);
  if (errno == ERANGE || text == end || *end != '\0' ||
      !std::isfinite(parsed)) {
    return false;
  }
  value = parsed;
  return true;
}

bool parse_parameters(const char* text, std::vector<fastf_t>& values,
                      std::size_t minimum)
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  double value = 0.0;
  while (input >> value) {
    if (!std::isfinite(value) || values.size() == RT_EDIT_MAXPARA) {
      return false;
    }
    values.push_back(value);
  }
  return input.eof() && values.size() >= minimum;
}

std::string vertex_label(std::size_t index, const point2d_t vertex,
                         double base2local)
{
  std::ostringstream text;
  text << index << "    " << std::setprecision(8) << vertex[0] * base2local
       << "    " << vertex[1] * base2local;
  return text.str();
}

std::string segment_label(std::size_t index, const void* segment)
{
  std::ostringstream text;
  text << index << "    ";
  if (!segment) {
    return text.str() + "Invalid";
  }
  const uint32_t magic = *static_cast<const uint32_t*>(segment);
  if (magic == CURVE_LSEG_MAGIC) {
    const auto* line = static_cast<const struct line_seg*>(segment);
    text << "Line    " << line->start << " -> " << line->end;
  } else if (magic == CURVE_CARC_MAGIC) {
    const auto* arc = static_cast<const struct carc_seg*>(segment);
    text << (arc->radius < 0.0 ? "Circle    " : "Arc    ") << arc->start
         << " -> " << arc->end << "    r=" << std::setprecision(8)
         << arc->radius;
  } else if (magic == CURVE_BEZIER_MAGIC) {
    const auto* bezier = static_cast<const struct bezier_seg*>(segment);
    text << "Bezier    degree=" << bezier->degree;
  } else if (magic == CURVE_NURB_MAGIC) {
    const auto* nurb = static_cast<const struct nurb_seg*>(segment);
    text << "NURB    order=" << nurb->order << " points=" << nurb->c_size;
  } else {
    text << "Unknown";
  }
  return text.str();
}

} // namespace

class Fl_Ged_Sketch_Editor::Impl
{
public:
  Impl(Fl_Ged_Sketch_Editor& owner, Fl_Cad_View::Backend backend)
    : owner_(owner)
  {
    owner_.begin();
    preview_ = new Fl_Cad_View(0, 0, owner_.w() - kPanelWidth, owner_.h(),
                               backend);
    preview_->draw_handler([this](struct bview* view) { draw(view); });
    preview_->error_handler(
      [this](const std::string& message) { report_error(message); });
    preview_->input_handler(
      [this](const Fl_Cad_View::InputEvent& event) {
        return handle_preview_input(event);
      });

    const int panel_x = owner_.w() - kPanelWidth + kMargin;
    heading_ = new Fl_Box(panel_x, kMargin, 174,
                          kRowHeight, "No sketch loaded");
    heading_->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
    pointer_mode_ = new Fl_Choice(panel_x + 180, kMargin, 144, kRowHeight);
    pointer_mode_->add("Navigate");
    pointer_mode_->add("Move Vertex");
    pointer_mode_->add("Add Vertex");
    pointer_mode_->value(0);
    pointer_mode_->tooltip(
      "Choose whether dragging navigates the preview or edits sketch vertices");
    vertices_ = new Fl_Hold_Browser(panel_x, kMargin + kRowHeight,
                                    kPanelWidth - 2 * kMargin, 180,
                                    "Vertices: index, U, V");
    vertices_->align(FL_ALIGN_TOP_LEFT);
    vertices_->callback(vertex_selected, this);

    const int vertex_y = kMargin + kRowHeight + 186;
    u_ = new Fl_Input(panel_x + 25, vertex_y, 120, kRowHeight, "U");
    v_ = new Fl_Input(panel_x + 174, vertex_y, 120, kRowHeight, "V");
    add_vertex_ = new Fl_Button(panel_x, vertex_y + kRowHeight + 4, 96,
                                kButtonHeight, "Add Vertex");
    move_vertex_ = new Fl_Button(panel_x + 103, vertex_y + kRowHeight + 4,
                                 106, kButtonHeight, "Move Vertex");
    delete_vertex_ = new Fl_Button(panel_x + 216,
                                   vertex_y + kRowHeight + 4, 108,
                                   kButtonHeight, "Delete Vertex");
    add_vertex_->callback(add_vertex_callback, this);
    move_vertex_->callback(move_vertex_callback, this);
    delete_vertex_->callback(delete_vertex_callback, this);

    const int segment_y = vertex_y + 2 * kRowHeight + kButtonHeight + 6;
    segments_ = new Fl_Hold_Browser(panel_x, segment_y,
                                    kPanelWidth - 2 * kMargin, 170,
                                    "Segments");
    segments_->align(FL_ALIGN_TOP_LEFT);
    start_ = new Fl_Input(panel_x + 48, segment_y + 176, 95, kRowHeight,
                          "Start");
    end_ = new Fl_Input(panel_x + 185, segment_y + 176, 95, kRowHeight,
                        "End");
    add_line_ = new Fl_Button(panel_x, segment_y + 207, 100, kButtonHeight,
                              "Add Line");
    add_arc_ = new Fl_Button(panel_x + 108, segment_y + 207, 100,
                             kButtonHeight, "Add Arc...");
    delete_segment_ = new Fl_Button(panel_x + 216, segment_y + 207, 108,
                                    kButtonHeight, "Delete Segment");
    add_bezier_ = new Fl_Button(panel_x, segment_y + 240, 156,
                                kButtonHeight, "Add Bezier...");
    add_nurb_ = new Fl_Button(panel_x + 164, segment_y + 240, 160,
                              kButtonHeight, "Add NURB...");
    add_line_->callback(add_line_callback, this);
    add_arc_->callback(add_arc_callback, this);
    add_bezier_->callback(add_bezier_callback, this);
    add_nurb_->callback(add_nurb_callback, this);
    delete_segment_->callback(delete_segment_callback, this);

    advanced_ = new Fl_Choice(panel_x + 82, segment_y + 273,
                              kPanelWidth - 2 * kMargin - 82, kRowHeight,
                              "Advanced");
    advanced_->add("Move vertex list...");
    advanced_->add("Constrain two vertices...");
    advanced_->add("Split selected segment...");
    advanced_->add("Reverse selected segment");
    advanced_->add("Toggle selected arc side");
    advanced_->add("Set selected arc radius...");
    advanced_->add("Set arc tangency...");
    advanced_->add("Edit selected NURB knots...");
    advanced_->add("Edit selected NURB weights...");
    advanced_->add("Set sketch plane...");
    advanced_->value(-1);
    advanced_->callback(advanced_callback, this);

    undo_ = new Fl_Button(panel_x, owner_.h() - kButtonHeight - kMargin, 88,
                          kButtonHeight, "Undo");
    save_ = new Fl_Button(panel_x + 98,
                          owner_.h() - kButtonHeight - kMargin, 88,
                          kButtonHeight, "Save");
    close_ = new Fl_Button(panel_x + 196,
                           owner_.h() - kButtonHeight - kMargin, 88,
                           kButtonHeight, "Close");
    undo_->callback(undo_callback, this);
    save_->callback(save_callback, this);
    close_->callback(close_callback, this);
    owner_.end();
    owner_.resizable(preview_);
    owner_.callback(close_callback, this);
    set_enabled(false);
  }

  ~Impl()
  {
    if (edit_) {
      rt_edit_destroy(edit_);
    }
  }

  bool open(struct ged* gedp, const std::string& object_name)
  {
    if (dirty_) {
      return fail("Save or discard the current sketch changes before opening another object");
    }
    reset();
    if (!gedp || !gedp->dbip || object_name.empty()) {
      return fail("Open a database and select a sketch first");
    }
    struct directory* candidate =
      db_lookup(gedp->dbip, object_name.c_str(), LOOKUP_QUIET);
    if (!candidate || candidate->d_minor_type != ID_SKETCH) {
      return fail("The selected object is not a sketch");
    }
    struct db_full_path path;
    db_full_path_init(&path);
    db_add_node_to_full_path(&path, candidate);
    struct rt_edit* candidate_edit =
      rt_edit_create(&path, gedp->dbip, &tolerance_, preview_->view());
    db_free_full_path(&path);
    if (!candidate_edit || candidate_edit->es_int.idb_type != ID_SKETCH) {
      if (candidate_edit) {
        rt_edit_destroy(candidate_edit);
      }
      return fail("BRL-CAD could not initialize the sketch editor");
    }
    auto* value = static_cast<struct rt_sketch_internal*>(
      candidate_edit->es_int.idb_ptr);
    if (value->vert_count > kMaxVertices || value->curve.count > kMaxSegments) {
      rt_edit_destroy(candidate_edit);
      return fail("The sketch exceeds the editor's vertex or segment limit");
    }

    ged_ = gedp;
    directory_ = candidate;
    edit_ = candidate_edit;
    object_name_ = object_name;
    saved_address_ = directory_->d_addr;
    saved_length_ = directory_->d_len;
    edit_->base2local = ged_->dbip->dbi_base2local;
    edit_->local2base = ged_->dbip->dbi_local2base;
    edit_->mv_context = 1;
    rt_edit_checkpoint(edit_);
    populate();
    fit_preview();
    heading_->copy_label(("Sketch: " + object_name_).c_str());
    owner_.copy_label(("Sketch Editor - " + object_name_).c_str());
    set_enabled(true);
    preview_->invalidate();
    last_error_.clear();
    return true;
  }

  bool save()
  {
    if (!edit_ || !ged_ || !ged_->dbip || !directory_) {
      return fail("No sketch is open for saving");
    }
    if (!dirty_) {
      return true;
    }
    if (ged_->dbip->dbi_read_only || db_version(ged_->dbip) < 5) {
      return fail("The database is read-only or does not support editing");
    }
    struct directory* current =
      db_lookup(ged_->dbip, object_name_.c_str(), LOOKUP_QUIET);
    if (current != directory_ || current->d_addr != saved_address_ ||
        current->d_len != saved_length_ || current->d_minor_type != ID_SKETCH) {
      return fail("The sketch changed outside this editor; reopen it before saving");
    }
    if (rt_db_put_internal(directory_, ged_->dbip, &edit_->es_int) < 0) {
      return fail("BRL-CAD could not write the edited sketch");
    }
    struct ged* context = ged_;
    const std::string saved_name = object_name_;
    dirty_ = false;
    checkpoint_dirty_ = false;
    if (!open(context, saved_name)) {
      return fail("Changes were saved, but the sketch edit session could not be rebuilt");
    }
    if (changed_handler_) {
      changed_handler_();
    }
    last_error_.clear();
    return true;
  }

  bool undo()
  {
    if (!edit_ || rt_edit_revert(edit_) != BRLCAD_OK) {
      return fail("There is no sketch checkpoint to restore");
    }
    dirty_ = checkpoint_dirty_;
    populate();
    preview_->invalidate();
    update_buttons();
    last_error_.clear();
    return true;
  }

  void resize(int width, int height)
  {
    const int panel_width = std::min(kPanelWidth, std::max(1, width / 2));
    const int panel_x = width - panel_width + kMargin;
    preview_->resize(0, 0, std::max(1, width - panel_width), height);
    const int content_width = panel_width - 2 * kMargin;
    heading_->resize(panel_x, kMargin, std::max(1, content_width - 150),
                     kRowHeight);
    pointer_mode_->resize(panel_x + std::max(1, content_width - 144), kMargin,
                          144, kRowHeight);
    vertices_->resize(panel_x, kMargin + kRowHeight, content_width, 180);
    const int vertex_y = kMargin + kRowHeight + 186;
    u_->resize(panel_x + 25, vertex_y, u_->w(), u_->h());
    v_->resize(panel_x + 174, vertex_y, v_->w(), v_->h());
    add_vertex_->position(panel_x, vertex_y + kRowHeight + 4);
    move_vertex_->position(panel_x + 103, vertex_y + kRowHeight + 4);
    delete_vertex_->position(panel_x + 216, vertex_y + kRowHeight + 4);
    const int segment_y = vertex_y + 2 * kRowHeight + kButtonHeight + 6;
    const int footer_y = height - kButtonHeight - kMargin;
    const int list_height = std::max(
      60, footer_y - segment_y - 3 * kRowHeight - 3 * kButtonHeight - 26);
    segments_->resize(panel_x, segment_y, content_width, list_height);
    start_->resize(panel_x + 48, segment_y + list_height + 6, start_->w(),
                   start_->h());
    end_->resize(panel_x + 185, segment_y + list_height + 6, end_->w(),
                 end_->h());
    add_line_->position(panel_x, segment_y + list_height + kRowHeight + 10);
    add_arc_->position(panel_x + 108,
                       segment_y + list_height + kRowHeight + 10);
    delete_segment_->position(panel_x + 216,
                              segment_y + list_height + kRowHeight + 10);
    add_bezier_->position(panel_x,
                          segment_y + list_height + kRowHeight +
                            kButtonHeight + 14);
    add_nurb_->position(panel_x + 164,
                        segment_y + list_height + kRowHeight +
                          kButtonHeight + 14);
    advanced_->resize(panel_x + 82,
                      segment_y + list_height + kRowHeight +
                        2 * kButtonHeight + 18,
                      std::max(1, content_width - 82), kRowHeight);
    undo_->position(panel_x, footer_y);
    save_->position(panel_x + 98, footer_y);
    close_->position(panel_x + 196, footer_y);
  }

  bool dirty() const noexcept { return dirty_; }
  const std::string& object_name() const noexcept { return object_name_; }
  const std::string& last_error() const noexcept { return last_error_; }

  bool move_vertex(std::size_t index, double u, double v)
  {
    const auto* value = sketch();
    if (!value || index >= value->vert_count || !std::isfinite(u) ||
        !std::isfinite(v)) {
      return fail("A valid vertex and finite U and V coordinates are required");
    }
    return pick_and_move_vertex(index, u, v);
  }

  bool move_vertices(const std::vector<std::size_t>& indices, double delta_u,
                     double delta_v)
  {
    const auto* value = sketch();
    if (!value || indices.empty() || indices.size() > RT_EDIT_MAXPARA - 2 ||
        !std::isfinite(delta_u) || !std::isfinite(delta_v)) {
      return fail("A bounded vertex list and finite movement are required");
    }
    std::vector<fastf_t> parameters{delta_u, delta_v};
    for (const auto index : indices) {
      if (index >= value->vert_count) {
        return fail("A vertex-list index is outside the sketch");
      }
      parameters.push_back(static_cast<fastf_t>(index));
    }
    return run_operation(ECMD_SKETCH_MOVE_VERTEX_LIST, std::move(parameters));
  }

  bool constrain_vertices(VertexConstraint constraint, std::size_t fixed,
                          std::size_t moving, double distance)
  {
    const auto* value = sketch();
    if (!value || fixed >= value->vert_count || moving >= value->vert_count ||
        fixed == moving) {
      return fail("Two distinct valid vertices are required for a constraint");
    }
    const double scale = edit_->base2local;
    const double fixed_u = value->verts[fixed][0] * scale;
    const double fixed_v = value->verts[fixed][1] * scale;
    double moving_u = value->verts[moving][0] * scale;
    double moving_v = value->verts[moving][1] * scale;
    switch (constraint) {
      case VertexConstraint::Coincident:
        moving_u = fixed_u;
        moving_v = fixed_v;
        break;
      case VertexConstraint::Horizontal: moving_v = fixed_v; break;
      case VertexConstraint::Vertical: moving_u = fixed_u; break;
      case VertexConstraint::Distance: {
        if (!std::isfinite(distance) || distance < 0.0) {
          return fail("Constraint distance must be finite and nonnegative");
        }
        double delta_u = moving_u - fixed_u;
        double delta_v = moving_v - fixed_v;
        const double length = std::hypot(delta_u, delta_v);
        if (length <= SMALL_FASTF) {
          delta_u = 1.0;
          delta_v = 0.0;
        } else {
          delta_u /= length;
          delta_v /= length;
        }
        moving_u = fixed_u + distance * delta_u;
        moving_v = fixed_v + distance * delta_v;
        break;
      }
    }
    return pick_and_move_vertex(moving, moving_u, moving_v);
  }

  bool split_segment(std::size_t index, double parameter)
  {
    const auto* value = sketch();
    if (!value || index >= value->curve.count || !std::isfinite(parameter) ||
        parameter <= 0.0 || parameter >= 1.0) {
      return fail("A valid segment and split parameter between zero and one are required");
    }
    return run_operation(ECMD_SKETCH_SPLIT_SEGMENT,
                         {static_cast<fastf_t>(index), parameter});
  }

  bool reverse_segment(std::size_t index)
  {
    return require_segment(index) &&
           run_operation(ECMD_SKETCH_TOGGLE_SEGMENT_REVERSE,
                         {static_cast<fastf_t>(index)});
  }

  bool set_arc_radius(std::size_t index, double radius)
  {
    if (!require_segment_type(index, CURVE_CARC_MAGIC, "arc")) {
      return false;
    }
    if (!std::isfinite(radius) || std::abs(radius) <= SMALL_FASTF) {
      return fail("Arc radius must be finite and nonzero");
    }
    return run_operation(ECMD_SKETCH_SET_ARC_RADIUS,
                         {static_cast<fastf_t>(index), radius});
  }

  bool toggle_arc_side(std::size_t index)
  {
    return require_segment_type(index, CURVE_CARC_MAGIC, "arc") &&
           run_operation(ECMD_SKETCH_TOGGLE_ARC_ORIENT,
                         {static_cast<fastf_t>(index)});
  }

  bool set_tangency(std::size_t arc, std::size_t adjacent,
                    double angle_radians)
  {
    if (!require_segment_type(arc, CURVE_CARC_MAGIC, "arc") ||
        !require_segment(adjacent) || arc == adjacent ||
        !std::isfinite(angle_radians)) {
      return fail("Tangency requires a valid arc, a different adjacent segment, and a finite angle");
    }
    return run_operation(ECMD_SKETCH_SET_TANGENCY,
                         {static_cast<fastf_t>(arc),
                          static_cast<fastf_t>(adjacent), angle_radians});
  }

  bool set_nurb_knots(std::size_t index, const std::vector<double>& knots)
  {
    if (!require_segment_type(index, CURVE_NURB_MAGIC, "NURB")) {
      return false;
    }
    const auto* nurb = static_cast<const struct nurb_seg*>(
      sketch()->curve.segment[index]);
    if (knots.size() != static_cast<std::size_t>(nurb->order + nurb->c_size) ||
        knots.size() + 1 > RT_EDIT_MAXPARA ||
        !std::is_sorted(knots.begin(), knots.end()) ||
        !std::all_of(knots.begin(), knots.end(),
                     [](double value) { return std::isfinite(value); })) {
      return fail("NURB knots must be finite, nondecreasing, and match order plus control-point count");
    }
    if (!pick_segment(index)) {
      return false;
    }
    std::vector<fastf_t> parameters{static_cast<fastf_t>(knots.size())};
    parameters.insert(parameters.end(), knots.begin(), knots.end());
    return run_operation(ECMD_SKETCH_NURB_EDIT_KV, std::move(parameters));
  }

  bool set_nurb_weights(std::size_t index,
                        const std::vector<double>& weights)
  {
    if (!require_segment_type(index, CURVE_NURB_MAGIC, "NURB")) {
      return false;
    }
    const auto* nurb = static_cast<const struct nurb_seg*>(
      sketch()->curve.segment[index]);
    if (weights.size() != static_cast<std::size_t>(nurb->c_size) ||
        weights.size() + 2 > RT_EDIT_MAXPARA ||
        !std::all_of(weights.begin(), weights.end(), [](double value) {
          return std::isfinite(value) && value > 0.0;
        })) {
      return fail("NURB weights must be positive finite values matching the control-point count");
    }
    std::vector<fastf_t> parameters{static_cast<fastf_t>(index),
                                    static_cast<fastf_t>(weights.size())};
    parameters.insert(parameters.end(), weights.begin(), weights.end());
    return run_operation(ECMD_SKETCH_NURB_EDIT_WEIGHTS,
                         std::move(parameters));
  }
  void changed_handler(ChangedHandler handler)
  {
    changed_handler_ = std::move(handler);
  }
  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

private:
  struct rt_sketch_internal* sketch() const
  {
    return edit_ ? static_cast<struct rt_sketch_internal*>(edit_->es_int.idb_ptr)
                 : nullptr;
  }

  void reset()
  {
    if (edit_) {
      rt_edit_destroy(edit_);
      edit_ = nullptr;
    }
    ged_ = nullptr;
    directory_ = nullptr;
    object_name_.clear();
    dirty_ = false;
    checkpoint_dirty_ = false;
    direct_dragging_ = false;
    direct_vertex_ = -1;
    pointer_mode_->value(0);
    vertices_->clear();
    segments_->clear();
    heading_->copy_label("No sketch loaded");
    set_enabled(false);
  }

  void populate()
  {
    auto* value = sketch();
    vertices_->clear();
    segments_->clear();
    if (!value) {
      return;
    }
    for (std::size_t index = 0; index < value->vert_count; ++index) {
      vertices_->add(
        vertex_label(index, value->verts[index], edit_->base2local).c_str());
    }
    for (std::size_t index = 0; index < value->curve.count; ++index) {
      segments_->add(segment_label(index, value->curve.segment[index]).c_str());
    }
  }

  bool run_operation(int command, std::vector<fastf_t> parameters)
  {
    if (!edit_ || parameters.size() > RT_EDIT_MAXPARA ||
        !EDOBJ[ID_SKETCH].ft_set_edit_mode) {
      return fail("The sketch operation is unavailable");
    }
    if (rt_edit_checkpoint(edit_) != BRLCAD_OK) {
      return fail("Unable to create an undo checkpoint");
    }
    checkpoint_dirty_ = dirty_;
    if (edit_->log_str) {
      bu_vls_trunc(edit_->log_str, 0);
    }
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, command);
    edit_->e_inpara = static_cast<int>(parameters.size());
    std::copy(parameters.begin(), parameters.end(), edit_->e_para);
    if (rt_edit_process(edit_) != BRLCAD_OK) {
      const std::string error = edit_log("BRL-CAD rejected the sketch edit");
      rt_edit_revert(edit_);
      return fail(error);
    }
    dirty_ = true;
    populate();
    preview_->invalidate();
    update_buttons();
    last_error_.clear();
    return true;
  }

  bool move_vertex()
  {
    const int selected = vertices_->value() - 1;
    double u = 0.0;
    double v = 0.0;
    if (selected < 0) {
      return fail("Select a vertex to move");
    }
    if (!parse_number(u_->value(), u) || !parse_number(v_->value(), v)) {
      return fail("Enter finite U and V coordinates");
    }
    return pick_and_move_vertex(static_cast<std::size_t>(selected), u, v);
  }

  bool pick_and_move_vertex(std::size_t selected, double u, double v)
  {
    if (rt_edit_checkpoint(edit_) != BRLCAD_OK) {
      return fail("Unable to create an undo checkpoint");
    }
    checkpoint_dirty_ = dirty_;
    if (edit_->log_str) {
      bu_vls_trunc(edit_->log_str, 0);
    }
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, ECMD_SKETCH_PICK_VERTEX);
    edit_->e_para[0] = static_cast<fastf_t>(selected);
    edit_->e_inpara = 1;
    if (rt_edit_process(edit_) != BRLCAD_OK) {
      return fail(edit_log("Unable to select the sketch vertex"));
    }
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, ECMD_SKETCH_MOVE_VERTEX);
    edit_->e_para[0] = u;
    edit_->e_para[1] = v;
    edit_->e_inpara = 2;
    if (rt_edit_process(edit_) != BRLCAD_OK) {
      const std::string error = edit_log("Unable to move the sketch vertex");
      rt_edit_revert(edit_);
      return fail(error);
    }
    dirty_ = true;
    populate();
    vertices_->value(static_cast<int>(selected) + 1);
    preview_->invalidate();
    update_buttons();
    return true;
  }

  bool add_vertex()
  {
    double u = 0.0;
    double v = 0.0;
    if (!parse_number(u_->value(), u) || !parse_number(v_->value(), v)) {
      return fail("Enter finite U and V coordinates");
    }
    return run_operation(ECMD_SKETCH_ADD_VERTEX, {u, v});
  }

  bool delete_vertex()
  {
    const int selected = vertices_->value() - 1;
    return selected < 0
             ? fail("Select a vertex to delete")
             : run_operation(ECMD_SKETCH_DELETE_VERTEX,
                             {static_cast<fastf_t>(selected)});
  }

  bool add_line()
  {
    double start = 0.0;
    double end = 0.0;
    if (!parse_number(start_->value(), start) ||
        !parse_number(end_->value(), end) || !EQUAL(std::trunc(start), start) ||
        !EQUAL(std::trunc(end), end)) {
      return fail("Line endpoints must be vertex indices");
    }
    return run_operation(ECMD_SKETCH_APPEND_LINE, {start, end});
  }

  bool add_arc()
  {
    const std::string defaults =
      std::string(start_->value() ? start_->value() : "0") + " " +
      (end_->value() ? end_->value() : "0") + " 1 1 0";
    const char* entered = fl_input(
      "start end radius center-is-left clockwise (0 or 1):",
      defaults.c_str());
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> parameters;
    if (!parse_parameters(entered, parameters, 5) || parameters.size() != 5) {
      return fail("An arc needs five finite values: start end radius side orientation");
    }
    return run_operation(ECMD_SKETCH_APPEND_ARC, std::move(parameters));
  }

  bool add_bezier()
  {
    const char* entered = fl_input(
      "Bezier control-point vertex indices (at least two):", "0 1");
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> parameters;
    if (!parse_parameters(entered, parameters, 2)) {
      return fail("Enter at least two valid control-point vertex indices");
    }
    return run_operation(ECMD_SKETCH_APPEND_BEZIER, std::move(parameters));
  }

  bool add_nurb()
  {
    const char* entered = fl_input(
      "NURB order followed by control-point vertex indices:", "2 0 1");
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> parameters;
    if (!parse_parameters(entered, parameters, 3)) {
      return fail("Enter an order and enough valid control-point indices");
    }
    return run_operation(ECMD_SKETCH_APPEND_NURB, std::move(parameters));
  }

  bool delete_segment()
  {
    const int selected = segments_->value() - 1;
    return selected < 0
             ? fail("Select a segment to delete")
             : run_operation(ECMD_SKETCH_DELETE_SEGMENT,
                             {static_cast<fastf_t>(selected)});
  }

  bool require_segment(std::size_t index)
  {
    const auto* value = sketch();
    return value && index < value->curve.count
             ? true
             : fail("The segment index is outside the sketch");
  }

  bool require_segment_type(std::size_t index, uint32_t magic,
                            const char* label)
  {
    if (!require_segment(index)) {
      return false;
    }
    const void* segment = sketch()->curve.segment[index];
    return segment && *static_cast<const uint32_t*>(segment) == magic
             ? true
             : fail("The selected segment is not a " + std::string(label));
  }

  bool pick_segment(std::size_t index)
  {
    if (!require_segment(index)) {
      return false;
    }
    if (edit_->log_str) {
      bu_vls_trunc(edit_->log_str, 0);
    }
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, ECMD_SKETCH_PICK_SEGMENT);
    edit_->e_para[0] = static_cast<fastf_t>(index);
    edit_->e_inpara = 1;
    return rt_edit_process(edit_) == BRLCAD_OK
             ? true
             : fail(edit_log("Unable to select the sketch segment"));
  }

  bool screen_to_uv(int x, int y, double& u, double& v) const
  {
    const auto* value = sketch();
    struct bview* view = preview_->view();
    if (!value || !view) {
      return false;
    }
    point_t model_point;
    if (bv_screen_pt(&model_point, x, y, view) != 0) {
      return false;
    }
    vect_t offset;
    VSUB2(offset, model_point, value->V);
    u = VDOT(offset, value->u_vec) * edit_->base2local;
    v = VDOT(offset, value->v_vec) * edit_->base2local;
    return std::isfinite(u) && std::isfinite(v);
  }

  int closest_vertex(int x, int y) const
  {
    const auto* value = sketch();
    struct bview* view = preview_->view();
    if (!value || !view || view->gv_width <= 0) {
      return -1;
    }
    fastf_t cursor_x = 0.0;
    fastf_t cursor_y = 0.0;
    if (bv_screen_to_view(view, &cursor_x, &cursor_y, x, y) != 0) {
      return -1;
    }
    const double radius = kPickRadiusPixels * 2.0 / view->gv_width;
    double best_distance = radius * radius;
    int best = -1;
    for (std::size_t index = 0; index < value->vert_count; ++index) {
      point_t model_point;
      VJOIN2(model_point, value->V, value->verts[index][0], value->u_vec,
             value->verts[index][1], value->v_vec);
      point_t view_point;
      MAT4X3PNT(view_point, view->gv_model2view, model_point);
      const double dx = view_point[X] - cursor_x;
      const double dy = view_point[Y] - cursor_y;
      const double distance = dx * dx + dy * dy;
      if (distance < best_distance) {
        best_distance = distance;
        best = static_cast<int>(index);
      }
    }
    return best;
  }

  bool begin_direct_move(int index)
  {
    if (index < 0 || rt_edit_checkpoint(edit_) != BRLCAD_OK) {
      return fail(index < 0 ? "No sketch vertex is near the pointer"
                            : "Unable to create an undo checkpoint");
    }
    checkpoint_dirty_ = dirty_;
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, ECMD_SKETCH_PICK_VERTEX);
    edit_->e_para[0] = index;
    edit_->e_inpara = 1;
    if (rt_edit_process(edit_) != BRLCAD_OK) {
      rt_edit_revert(edit_);
      return fail(edit_log("Unable to select the sketch vertex"));
    }
    direct_vertex_ = index;
    direct_dragging_ = true;
    vertices_->value(index + 1);
    vertex_selected(nullptr, this);
    return true;
  }

  bool continue_direct_move(int x, int y)
  {
    double u = 0.0;
    double v = 0.0;
    if (!direct_dragging_ || !screen_to_uv(x, y, u, v)) {
      return false;
    }
    EDOBJ[ID_SKETCH].ft_set_edit_mode(edit_, ECMD_SKETCH_MOVE_VERTEX);
    edit_->e_para[0] = u;
    edit_->e_para[1] = v;
    edit_->e_inpara = 2;
    if (rt_edit_process(edit_) != BRLCAD_OK) {
      const std::string error = edit_log("Unable to move the sketch vertex");
      rt_edit_revert(edit_);
      direct_dragging_ = false;
      return fail(error);
    }
    dirty_ = true;
    populate();
    vertices_->value(direct_vertex_ + 1);
    vertex_selected(nullptr, this);
    preview_->invalidate();
    update_buttons();
    return true;
  }

  bool handle_preview_input(const Fl_Cad_View::InputEvent& event)
  {
    if (!edit_ || !pointer_mode_ || pointer_mode_->value() <= 0) {
      return false;
    }
    if (event.button != FL_LEFT_MOUSE &&
        (event.buttons & FL_BUTTON1) == 0 && event.type != FL_RELEASE) {
      return true;
    }
    if (pointer_mode_->value() == 2 && event.type == FL_PUSH) {
      double u = 0.0;
      double v = 0.0;
      return screen_to_uv(event.x, event.y, u, v) &&
             run_operation(ECMD_SKETCH_ADD_VERTEX, {u, v});
    }
    if (pointer_mode_->value() == 1) {
      if (event.type == FL_PUSH) {
        begin_direct_move(closest_vertex(event.x, event.y));
        return true;
      }
      if (event.type == FL_DRAG) {
        continue_direct_move(event.x, event.y);
        return true;
      }
      if (event.type == FL_RELEASE) {
        direct_dragging_ = false;
        direct_vertex_ = -1;
        return true;
      }
    }
    return true;
  }

  bool advanced_move_vertices()
  {
    const char* entered = fl_input(
      "U-delta V-delta followed by vertex indices (maximum 18):", "0 0 0");
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> parameters;
    if (!parse_parameters(entered, parameters, 3)) {
      return fail("Enter two finite deltas followed by at least one vertex index");
    }
    std::vector<std::size_t> indices;
    indices.reserve(parameters.size() - 2);
    for (std::size_t offset = 2; offset < parameters.size(); ++offset) {
      if (parameters[offset] < 0.0 ||
          !EQUAL(std::trunc(parameters[offset]), parameters[offset])) {
        return fail("Vertex-list indices must be nonnegative integers");
      }
      indices.push_back(static_cast<std::size_t>(parameters[offset]));
    }
    return move_vertices(indices, parameters[0], parameters[1]);
  }

  bool advanced_constraint()
  {
    const char* entered = fl_input(
      "Constraint c|h|v|d, fixed vertex, moving vertex, optional distance:",
      "h 0 1");
    if (!entered) {
      return true;
    }
    std::istringstream input(entered);
    char type = '\0';
    long long fixed = -1;
    long long moving = -1;
    double distance = 0.0;
    if (!(input >> type >> fixed >> moving) || fixed < 0 || moving < 0) {
      return fail("Enter a constraint type and two nonnegative vertex indices");
    }
    VertexConstraint constraint = VertexConstraint::Coincident;
    switch (type) {
      case 'c': constraint = VertexConstraint::Coincident; break;
      case 'h': constraint = VertexConstraint::Horizontal; break;
      case 'v': constraint = VertexConstraint::Vertical; break;
      case 'd':
        constraint = VertexConstraint::Distance;
        if (!(input >> distance)) {
          return fail("A distance constraint also requires a distance");
        }
        break;
      default: return fail("Constraint type must be c, h, v, or d");
    }
    input >> std::ws;
    if (!input.eof()) {
      return fail("The constraint contains unexpected values");
    }
    return constrain_vertices(constraint, static_cast<std::size_t>(fixed),
                              static_cast<std::size_t>(moving), distance);
  }

  bool advanced_split()
  {
    const int selected = segments_->value() - 1;
    if (selected < 0) {
      return fail("Select a segment to split");
    }
    const char* entered = fl_input("Split parameter between 0 and 1:", "0.5");
    if (!entered) {
      return true;
    }
    double parameter = 0.0;
    if (!parse_number(entered, parameter)) {
      return fail("Enter a finite split parameter between 0 and 1");
    }
    return split_segment(static_cast<std::size_t>(selected), parameter);
  }

  bool advanced_reverse()
  {
    const int selected = segments_->value() - 1;
    return selected >= 0
             ? reverse_segment(static_cast<std::size_t>(selected))
             : fail("Select a segment to reverse");
  }

  bool advanced_arc_side()
  {
    const int selected = segments_->value() - 1;
    return selected >= 0
             ? toggle_arc_side(static_cast<std::size_t>(selected))
             : fail("Select an arc first");
  }

  bool advanced_arc_radius()
  {
    const int selected = segments_->value() - 1;
    if (selected < 0) {
      return fail("Select an arc first");
    }
    const char* entered = fl_input("New arc radius in local units:", "1");
    if (!entered) {
      return true;
    }
    double radius = 0.0;
    if (!parse_number(entered, radius)) {
      return fail("Enter a finite nonzero arc radius");
    }
    return set_arc_radius(static_cast<std::size_t>(selected), radius);
  }

  bool advanced_tangency()
  {
    const int selected = segments_->value() - 1;
    if (selected < 0) {
      return fail("Select the arc to constrain first");
    }
    const char* entered =
      fl_input("Adjacent segment index and angle in radians:", "0 0");
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> values;
    if (!parse_parameters(entered, values, 2) || values.size() != 2 ||
        values[0] < 0.0 || !EQUAL(std::trunc(values[0]), values[0])) {
      return fail("Enter an adjacent segment index and finite angle");
    }
    return set_tangency(static_cast<std::size_t>(selected),
                        static_cast<std::size_t>(values[0]), values[1]);
  }

  bool advanced_nurb_knots()
  {
    const int selected = segments_->value() - 1;
    if (selected < 0 ||
        !require_segment_type(static_cast<std::size_t>(selected),
                              CURVE_NURB_MAGIC, "NURB")) {
      return selected < 0 ? fail("Select a NURB segment first") : false;
    }
    const auto* nurb = static_cast<const struct nurb_seg*>(
      sketch()->curve.segment[selected]);
    std::ostringstream defaults;
    for (int index = 0; index < nurb->k.k_size; ++index) {
      defaults << (index ? " " : "") << nurb->k.knots[index];
    }
    const char* entered = fl_input("Replacement nondecreasing knot values:",
                                   defaults.str().c_str());
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> values;
    if (!parse_parameters(entered, values,
                          static_cast<std::size_t>(nurb->k.k_size)) ||
        values.size() != static_cast<std::size_t>(nurb->k.k_size)) {
      return fail("Enter exactly order plus control-point-count knot values");
    }
    return set_nurb_knots(static_cast<std::size_t>(selected),
                          std::vector<double>(values.begin(), values.end()));
  }

  bool advanced_nurb_weights()
  {
    const int selected = segments_->value() - 1;
    if (selected < 0 ||
        !require_segment_type(static_cast<std::size_t>(selected),
                              CURVE_NURB_MAGIC, "NURB")) {
      return selected < 0 ? fail("Select a NURB segment first") : false;
    }
    const auto* nurb = static_cast<const struct nurb_seg*>(
      sketch()->curve.segment[selected]);
    std::ostringstream defaults;
    for (int index = 0; index < nurb->c_size; ++index) {
      defaults << (index ? " " : "")
               << (nurb->weights ? nurb->weights[index] : 1.0);
    }
    const char* entered = fl_input("Positive NURB control-point weights:",
                                   defaults.str().c_str());
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> values;
    if (!parse_parameters(entered, values,
                          static_cast<std::size_t>(nurb->c_size)) ||
        values.size() != static_cast<std::size_t>(nurb->c_size)) {
      return fail("Enter exactly one positive weight per NURB control point");
    }
    return set_nurb_weights(static_cast<std::size_t>(selected),
                            std::vector<double>(values.begin(), values.end()));
  }

  bool advanced_plane()
  {
    const auto* value = sketch();
    if (!value) {
      return fail("Open a sketch before changing its plane");
    }
    std::ostringstream defaults;
    defaults << value->V[0] * edit_->base2local << " "
             << value->V[1] * edit_->base2local << " "
             << value->V[2] * edit_->base2local << " " << value->u_vec[0]
             << " " << value->u_vec[1] << " " << value->u_vec[2] << " "
             << value->v_vec[0] << " " << value->v_vec[1] << " "
             << value->v_vec[2];
    const char* entered = fl_input(
      "Origin X Y Z, U direction X Y Z, V direction X Y Z:",
      defaults.str().c_str());
    if (!entered) {
      return true;
    }
    std::vector<fastf_t> parameters;
    if (!parse_parameters(entered, parameters, 9) || parameters.size() != 9) {
      return fail("The sketch plane requires nine finite values");
    }
    return run_operation(ECMD_SKETCH_SET_PLANE, std::move(parameters));
  }

  void fit_preview()
  {
    if (!edit_ || !preview_->view()) {
      return;
    }
    point_t minimum;
    point_t maximum;
    if (!OBJ[ID_SKETCH].ft_bbox ||
        OBJ[ID_SKETCH].ft_bbox(&edit_->es_int, &minimum, &maximum,
                               &tolerance_) < 0) {
      return;
    }
    point_t center;
    VADD2SCALE(center, minimum, maximum, 0.5);
    MAT_DELTAS_VEC_NEG(preview_->view()->gv_center, center);
    preview_->view()->gv_scale =
      std::max(1.0, static_cast<double>(DIST_PNT_PNT(minimum, maximum)) * 0.65);
    preview_->view()->gv_size = 2.0 * preview_->view()->gv_scale;
    preview_->view()->gv_isize = 1.0 / preview_->view()->gv_size;
    bv_update(preview_->view());
  }

  void draw(struct bview* view)
  {
    if (!edit_ || !view || !view->dmp) {
      return;
    }
    struct bu_list vlist;
    BU_LIST_INIT(&vlist);
    struct bg_tess_tol tessellation;
    tessellation.magic = BG_TESS_TOL_MAGIC;
    tessellation.abs = 0.0;
    tessellation.rel = 0.01;
    tessellation.norm = 0.0;
    if (OBJ[ID_SKETCH].ft_plot(
          &vlist, &edit_->es_int, &tessellation, &tolerance_, view) == 0) {
      auto* display = static_cast<struct dm*>(view->dmp);
      dm_set_fg(display, kPreviewRed, kPreviewGreen, kPreviewBlue, 1, 1.0);
      dm_draw_vlist(display, reinterpret_cast<struct bv_vlist*>(&vlist));
    }
    bv_vlist_cleanup(&vlist);
  }

  std::string edit_log(const char* fallback) const
  {
    if (edit_ && edit_->log_str && bu_vls_strlen(edit_->log_str) > 0) {
      const std::size_t length = bu_vls_strlen(edit_->log_str);
      const std::size_t accepted = std::min(length, kMaxLogBytes);
      std::string result(bu_vls_cstr(edit_->log_str), accepted);
      if (accepted != length) {
        result += "\n[edit diagnostic truncated]";
      }
      return result;
    }
    return fallback;
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* widgets[] = {vertices_, u_, v_, add_vertex_, move_vertex_,
                            delete_vertex_, segments_, start_, end_, add_line_,
                            add_arc_, add_bezier_, add_nurb_, delete_segment_,
                            advanced_, pointer_mode_, undo_, save_};
    for (auto* widget : widgets) {
      if (enabled) {
        widget->activate();
      } else {
        widget->deactivate();
      }
    }
    update_buttons();
  }

  void update_buttons()
  {
    if (edit_ && dirty_) {
      save_->activate();
    } else {
      save_->deactivate();
    }
  }

  bool fail(std::string message)
  {
    report_error(std::move(message));
    return false;
  }

  void report_error(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
  }

  bool request_close()
  {
    if (!dirty_) {
      owner_.hide();
      return true;
    }
    const int result = fl_choice("Save changes to %s before closing?", "Cancel",
                                 "Discard", "Save", object_name_.c_str());
    if (result == 0) {
      return false;
    }
    if (result == 2 && !save()) {
      fl_alert("%s", last_error_.c_str());
      return false;
    }
    if (result == 1) {
      dirty_ = false;
      checkpoint_dirty_ = false;
      update_buttons();
    }
    owner_.hide();
    return true;
  }

  static void vertex_selected(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    const int selected = self.vertices_->value() - 1;
    const auto* value = self.sketch();
    if (selected < 0 || !value ||
        static_cast<std::size_t>(selected) >= value->vert_count) {
      return;
    }
    const double scale = self.edit_->base2local;
    self.u_->value(std::to_string(value->verts[selected][0] * scale).c_str());
    self.v_->value(std::to_string(value->verts[selected][1] * scale).c_str());
  }

  template <bool (Impl::*Operation)()>
  static void operation_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!(self.*Operation)()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }

  static void add_vertex_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::add_vertex>(widget, data);
  }
  static void move_vertex_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::move_vertex>(widget, data);
  }
  static void delete_vertex_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::delete_vertex>(widget, data);
  }
  static void add_line_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::add_line>(widget, data);
  }
  static void add_arc_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::add_arc>(widget, data);
  }
  static void add_bezier_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::add_bezier>(widget, data);
  }
  static void add_nurb_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::add_nurb>(widget, data);
  }
  static void delete_segment_callback(Fl_Widget* widget, void* data)
  {
    operation_callback<&Impl::delete_segment>(widget, data);
  }
  static void advanced_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    using Operation = bool (Impl::*)();
    const Operation operations[] = {
      &Impl::advanced_move_vertices, &Impl::advanced_constraint,
      &Impl::advanced_split,         &Impl::advanced_reverse,
      &Impl::advanced_arc_side,      &Impl::advanced_arc_radius,
      &Impl::advanced_tangency,      &Impl::advanced_nurb_knots,
      &Impl::advanced_nurb_weights,  &Impl::advanced_plane};
    const int selected = self.advanced_->value();
    self.advanced_->value(-1);
    if (selected >= 0 &&
        static_cast<std::size_t>(selected) < std::size(operations) &&
        !(self.*operations[selected])()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }
  static void undo_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.undo()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }
  static void save_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.save()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }
  static void close_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->request_close();
  }

  Fl_Ged_Sketch_Editor& owner_;
  Fl_Cad_View* preview_ = nullptr;
  Fl_Box* heading_ = nullptr;
  Fl_Choice* pointer_mode_ = nullptr;
  Fl_Hold_Browser* vertices_ = nullptr;
  Fl_Input* u_ = nullptr;
  Fl_Input* v_ = nullptr;
  Fl_Button* add_vertex_ = nullptr;
  Fl_Button* move_vertex_ = nullptr;
  Fl_Button* delete_vertex_ = nullptr;
  Fl_Hold_Browser* segments_ = nullptr;
  Fl_Input* start_ = nullptr;
  Fl_Input* end_ = nullptr;
  Fl_Button* add_line_ = nullptr;
  Fl_Button* add_arc_ = nullptr;
  Fl_Button* add_bezier_ = nullptr;
  Fl_Button* add_nurb_ = nullptr;
  Fl_Button* delete_segment_ = nullptr;
  Fl_Choice* advanced_ = nullptr;
  Fl_Button* undo_ = nullptr;
  Fl_Button* save_ = nullptr;
  Fl_Button* close_ = nullptr;
  struct ged* ged_ = nullptr;
  struct directory* directory_ = nullptr;
  struct rt_edit* edit_ = nullptr;
  struct bn_tol tolerance_ = BN_TOL_INIT_TOL;
  b_off_t saved_address_ = 0;
  size_t saved_length_ = 0;
  std::string object_name_;
  std::string last_error_;
  bool dirty_ = false;
  bool checkpoint_dirty_ = false;
  bool direct_dragging_ = false;
  int direct_vertex_ = -1;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_Sketch_Editor::Fl_Ged_Sketch_Editor(Fl_Cad_View::Backend backend)
  : Fl_Double_Window(kInitialWidth, kInitialHeight, "Sketch Editor"),
    impl_(std::make_unique<Impl>(*this, backend))
{
}
Fl_Ged_Sketch_Editor::~Fl_Ged_Sketch_Editor() = default;
bool Fl_Ged_Sketch_Editor::open(struct ged* gedp,
                                const std::string& object_name)
{
  return impl_->open(gedp, object_name);
}
bool Fl_Ged_Sketch_Editor::save() { return impl_->save(); }
bool Fl_Ged_Sketch_Editor::undo() { return impl_->undo(); }
bool Fl_Ged_Sketch_Editor::move_vertex(std::size_t index, double u, double v)
{
  return impl_->move_vertex(index, u, v);
}
bool Fl_Ged_Sketch_Editor::move_vertices(
  const std::vector<std::size_t>& indices, double delta_u, double delta_v)
{
  return impl_->move_vertices(indices, delta_u, delta_v);
}
bool Fl_Ged_Sketch_Editor::constrain_vertices(VertexConstraint constraint,
                                              std::size_t fixed,
                                              std::size_t moving,
                                              double distance)
{
  return impl_->constrain_vertices(constraint, fixed, moving, distance);
}
bool Fl_Ged_Sketch_Editor::split_segment(std::size_t index, double parameter)
{
  return impl_->split_segment(index, parameter);
}
bool Fl_Ged_Sketch_Editor::reverse_segment(std::size_t index)
{
  return impl_->reverse_segment(index);
}
bool Fl_Ged_Sketch_Editor::set_arc_radius(std::size_t index, double radius)
{
  return impl_->set_arc_radius(index, radius);
}
bool Fl_Ged_Sketch_Editor::toggle_arc_side(std::size_t index)
{
  return impl_->toggle_arc_side(index);
}
bool Fl_Ged_Sketch_Editor::set_tangency(std::size_t arc,
                                        std::size_t adjacent,
                                        double angle_radians)
{
  return impl_->set_tangency(arc, adjacent, angle_radians);
}
bool Fl_Ged_Sketch_Editor::set_nurb_knots(
  std::size_t index, const std::vector<double>& knots)
{
  return impl_->set_nurb_knots(index, knots);
}
bool Fl_Ged_Sketch_Editor::set_nurb_weights(
  std::size_t index, const std::vector<double>& weights)
{
  return impl_->set_nurb_weights(index, weights);
}
bool Fl_Ged_Sketch_Editor::dirty() const noexcept { return impl_->dirty(); }
const std::string& Fl_Ged_Sketch_Editor::object_name() const noexcept
{
  return impl_->object_name();
}
const std::string& Fl_Ged_Sketch_Editor::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Sketch_Editor::changed_handler(ChangedHandler handler)
{
  impl_->changed_handler(std::move(handler));
}
void Fl_Ged_Sketch_Editor::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
}
void Fl_Ged_Sketch_Editor::resize(int x, int y, int width, int height)
{
  Fl_Double_Window::resize(x, y, width, height);
  if (impl_) {
    impl_->resize(width, height);
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
