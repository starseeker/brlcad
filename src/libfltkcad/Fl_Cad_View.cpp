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

#include "fltkcad/Fl_Cad_View.H"

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <FL/fl_draw.H>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

extern "C" {
#include "bn.h"
#include "bu/malloc.h"
#include "bu/vls.h"
#include "bv/util.h"
#include "dm.h"
#include "dm/extgl.h"
#include "dm/faceplate.h"
#include "ged/edit.h"
#include "ged/view.h"
#include "vmath.h"
}

namespace {

/* Substitute preview vectors only for the duration of a frame.  Keeping the
 * database scene intact lets reject, redraw, and view-specific LoD caches
 * retain their normal ownership and invalidation rules. */
class EditPreview
{
public:
  EditPreview(struct ged* context, struct bview* view)
  {
    if (!context || ged_edit_session_mode(context) == GED_EDIT_VIEW) return;
    for (int flags : {BV_DB_OBJS, BV_DB_OBJS | BV_LOCAL_OBJS}) {
      struct bu_ptbl* objects = bv_view_objs(view, flags);
      if (!objects) continue;
      for (size_t i = 0; i < BU_PTBL_LEN(objects); ++i)
        prepare(context, view, reinterpret_cast<struct bv_scene_obj*>(BU_PTBL_GET(objects, i)));
    }
  }

private:
  struct Vectors {
    explicit Vectors(struct bv_scene_obj* object) : object_(object), flags_(object->s_type_flags),
      update_(object->s_update_callback), length_(object->s_vlen)
    {
      BU_LIST_INIT(&saved_);
      BU_LIST_APPEND_LIST(&saved_, &object_->s_vlist);
      object_->s_type_flags &= ~(BV_MESH_LOD | BV_CSG_LOD);
      object_->s_update_callback = nullptr;
    }
    ~Vectors()
    {
      bv_vlist_cleanup(&object_->s_vlist);
      BU_LIST_APPEND_LIST(&object_->s_vlist, &saved_);
      object_->s_type_flags = flags_;
      object_->s_update_callback = update_;
      object_->s_vlen = length_;
    }
    struct bv_scene_obj* object_;
    struct bu_list saved_;
    unsigned long long flags_;
    int (*update_)(struct bv_scene_obj*, struct bview*, int);
    size_t length_;
  };

  void prepare(struct ged* context, struct bview* view, struct bv_scene_obj* object)
  {
    if (std::find(visited_.begin(), visited_.end(), object) != visited_.end()) return;
    visited_.push_back(object);
    for (size_t i = 0; i < BU_PTBL_LEN(&object->children); ++i)
      prepare(context, view, reinterpret_cast<struct bv_scene_obj*>(BU_PTBL_GET(&object->children, i)));
    if (!(object->s_type_flags & BV_DB_OBJS)) return;
    const struct db_full_path* path = static_cast<struct db_full_path*>(object->s_path);
    if (!path && object->s_u_data)
      path = &static_cast<struct ged_bv_data*>(object->s_u_data)->s_fullpath;
    if (!path) return;
    struct bu_list preview;
    BU_LIST_INIT(&preview);
    int ret = ged_edit_session_plot(context, path, view, &preview);
    if (ret == 1) {
      struct bv_scene_obj* rendered = bv_obj_for_view(object, view);
      if (!rendered) rendered = object;
      replacements_.emplace_back(std::make_unique<Vectors>(rendered));
      BU_LIST_APPEND_LIST(&rendered->s_vlist, &preview);
      rendered->s_vlen = bv_vlist_cmd_cnt(reinterpret_cast<struct bv_vlist*>(&rendered->s_vlist));
    }
    bv_vlist_cleanup(&preview);
  }
  std::vector<struct bv_scene_obj*> visited_;
  std::vector<std::unique_ptr<Vectors>> replacements_;
};

constexpr fastf_t kViewMinimumDepth = -100.0;
constexpr fastf_t kViewMaximumDepth = 100.0;
constexpr double kMouseZoomOut = 0.5;
constexpr double kMouseZoomIn = 2.0;
constexpr double kWheelZoomIn = 1.1;
constexpr double kWheelZoomOut = 0.9;
constexpr double kRotationDegreesPerPixel = 0.25;
constexpr double kScalePixels = 200.0;
constexpr unsigned char kErrorBackground = 32;
constexpr int kOSMesaMaxFramebufferDimension = 4096;
constexpr int kMaxGlErrorDrain = 16;

bool valid_dimensions(int width, int height, const Fl_Cad_View::Limits& limits)
{
  if (width <= 0 || height <= 0 || width > limits.max_framebuffer_width ||
      height > limits.max_framebuffer_height) {
    return false;
  }
  return static_cast<std::size_t>(width) <=
         limits.max_framebuffer_pixels / static_cast<std::size_t>(height);
}

bool valid_backend_dimensions(int width, int height,
                              Fl_Cad_View::Backend backend,
                              const Fl_Cad_View::Limits& limits)
{
  if (!valid_dimensions(width, height, limits)) {
    return false;
  }
  return backend != Fl_Cad_View::Backend::OSMesa ||
         (width <= kOSMesaMaxFramebufferDimension &&
          height <= kOSMesaMaxFramebufferDimension);
}

void flip_rows(unsigned char* pixels, int width, int height, int channels)
{
  const std::size_t row_bytes = static_cast<std::size_t>(width) * channels;
  std::vector<unsigned char> row(row_bytes);
  for (int y = 0; y < height / 2; ++y) {
    auto* top = pixels + static_cast<std::size_t>(y) * row_bytes;
    auto* bottom = pixels + static_cast<std::size_t>(height - y - 1) * row_bytes;
    std::memcpy(row.data(), top, row_bytes);
    std::memcpy(top, bottom, row_bytes);
    std::memcpy(bottom, row.data(), row_bytes);
  }
}

} // namespace

class Fl_Cad_View::Impl
{
public:
  class Canvas : public Fl_Widget
  {
  public:
    Canvas(Impl& owner, int x, int y, int width, int height)
      : Fl_Widget(x, y, width, height), owner_(owner)
    {
    }

    virtual Backend backend() const noexcept = 0;
    virtual int framebuffer_width() const = 0;
    virtual int framebuffer_height() const = 0;
    virtual bool context_ready() const { return true; }
    virtual void make_host_current() {}

  protected:
    Impl& owner_;
  };

  class SoftwareCanvas final : public Canvas
  {
  public:
    using Canvas::Canvas;

    Backend backend() const noexcept override { return Backend::OSMesa; }
    int framebuffer_width() const override { return w(); }
    int framebuffer_height() const override { return h(); }

    void draw() override
    {
      if (!owner_.render()) {
        owner_.draw_error(*this);
        return;
      }
      if (owner_.pixels_.empty()) {
        return;
      }
      fl_draw_image(owner_.pixels_.data(), x(), y(), owner_.pixel_width_,
                    owner_.pixel_height_, 4, owner_.pixel_width_ * 4);
    }

    int handle(int event) override
    {
      return owner_.handle_input(event) ? 1 : Canvas::handle(event);
    }
  };

  class GlCanvas final : public Fl_Gl_Window
  {
  public:
    GlCanvas(Impl& owner, int x, int y, int width, int height)
      : Fl_Gl_Window(x, y, width, height), owner_(owner)
    {
      mode(FL_RGB | FL_ALPHA | FL_DOUBLE | FL_DEPTH);
    }

    Backend backend() const noexcept { return Backend::SystemOpenGL; }
    int framebuffer_width() { return pixel_w(); }
    int framebuffer_height() { return pixel_h(); }
    bool ensure_context()
    {
      if (!shown()) {
        return false;
      }
      make_current();
      return context() != nullptr;
    }
    void make_host_current() { make_current(); }

    void draw() override
    {
      if (!owner_.render()) {
        if (!owner_.request_software_fallback()) {
          owner_.draw_error(*this);
        }
      }
    }

    int handle(int event) override
    {
      return owner_.handle_input(event) ? 1 : Fl_Gl_Window::handle(event);
    }

  private:
    Impl& owner_;
  };

  Impl(Fl_Cad_View& owner, Backend requested)
    : owner_(owner), preference_(requested)
  {
    BU_GET(local_view_, struct bview);
    bv_init(local_view_, nullptr);
    bu_vls_sprintf(&local_view_->gv_name, "fltk");
    active_view_ = local_view_;

    Backend selected = requested;
    if (selected == Backend::Auto) {
      const int mode = FL_RGB | FL_ALPHA | FL_DOUBLE | FL_DEPTH;
      selected = Fl_Gl_Window::can_do(mode) ? Backend::SystemOpenGL
                                            : Backend::OSMesa;
    }

    owner_.begin();
    if (selected == Backend::SystemOpenGL) {
      auto* gl_canvas = new GlCanvas(*this, owner_.x(), owner_.y(), owner_.w(),
                                     owner_.h());
      gl_canvas_ = gl_canvas;
      canvas_ = gl_canvas;
    } else {
      canvas_ = new SoftwareCanvas(*this, owner_.x(), owner_.y(), owner_.w(),
                                   owner_.h());
    }
    owner_.end();
    owner_.resizable(canvas_);
  }

  ~Impl()
  {
    if (fallback_pending_) {
      Fl::remove_timeout(software_fallback, this);
    }
    close_display_manager();
    if (local_view_) {
      bv_free(local_view_);
      BU_PUT(local_view_, struct bview);
    }
  }

  Backend actual_backend() const noexcept
  {
    return gl_canvas_ ? Backend::SystemOpenGL : Backend::OSMesa;
  }

  bool request_software_fallback()
  {
    if (preference_ != Backend::Auto || !gl_canvas_ || fallback_pending_) {
      return false;
    }
    fallback_pending_ = true;
    Fl::add_timeout(0.0, software_fallback, this);
    return true;
  }

  static void software_fallback(void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    self.fallback_pending_ = false;
    if (!self.gl_canvas_) {
      return;
    }

    self.close_display_manager();
    self.owner_.remove(self.gl_canvas_);
    delete self.gl_canvas_;
    self.gl_canvas_ = nullptr;

    self.owner_.begin();
    self.canvas_ = new SoftwareCanvas(
      self, self.owner_.x(), self.owner_.y(), self.owner_.w(), self.owner_.h());
    self.owner_.end();
    self.owner_.resizable(self.canvas_);
    self.pixels_.clear();
    self.pixel_width_ = 0;
    self.pixel_height_ = 0;
    self.dirty_ = true;
    self.canvas_->redraw();
  }

  bool initialize()
  {
    if (dmp_) {
      return true;
    }
    last_error_.clear();

    const int width = framebuffer_width();
    const int height = framebuffer_height();
    if (!valid_backend_dimensions(width, height, actual_backend(), limits_)) {
      report_error("Display framebuffer dimensions exceed configured limits");
      return false;
    }

    active_view_->gv_width = width;
    active_view_->gv_height = height;
    const char* attach = "attach";
    if (actual_backend() == Backend::SystemOpenGL) {
      if (!gl_canvas_ || !gl_canvas_->ensure_context()) {
        report_error("OpenGL display cannot initialize before its window is shown");
        return false;
      }
      extgl_context_ = {};
      extgl_context_.api_version = DM_EXTGL_API;
      extgl_context_.data = this;
      extgl_context_.make_current = extgl_make_current;
      extgl_context_.get_framebuffer_size = extgl_get_size;
      extgl_context_.get_display_image = extgl_capture;
      extgl_context_.draw_string = extgl_draw_string;
      extgl_context_.string_bbox = extgl_string_bbox;
      dmp_ = dm_open(&extgl_context_, nullptr, "extgl", 1, &attach);
    } else {
      dmp_ = dm_open(active_view_, nullptr, "swrast", 1, &attach);
    }
    if (!dmp_) {
      report_error(actual_backend() == Backend::SystemOpenGL
                     ? "Unable to open the libdm extgl backend"
                     : "Unable to open the libdm OSMesa swrast backend");
      return false;
    }

    dm_set_pathname(dmp_, actual_backend() == Backend::SystemOpenGL
                            ? "FLTKGL"
                            : "FLTKSW");
    fastf_t bounds[6] = {-1.0, 1.0, -1.0, 1.0,
                         kViewMinimumDepth, kViewMaximumDepth};
    dm_set_win_bounds(dmp_, bounds);
    dm_set_vp(dmp_, &active_view_->gv_scale);
    active_view_->dmp = dmp_;
    dm_set_faceplate_gui(dmp_, faceplate_gui_);
    dm_set_zbuffer(dmp_, 1);
    if (dm_set_bg(dmp_, background_first_.red, background_first_.green,
                  background_first_.blue, background_second_.red,
                  background_second_.green, background_second_.blue) !=
        BRLCAD_OK) {
      report_error("Unable to apply the display background color");
      close_display_manager();
      return false;
    }
    if (dm_configure_win(dmp_, 1) != BRLCAD_OK) {
      report_error("Unable to configure the libdm framebuffer");
      close_display_manager();
      return false;
    }
    dirty_ = true;
    return true;
  }

  void close_display_manager()
  {
    if (!dmp_) {
      return;
    }
    if (active_view_ && active_view_->dmp == dmp_) {
      active_view_->dmp = nullptr;
    }
    dm_close(dmp_);
    dmp_ = nullptr;
  }

  bool render()
  {
    if (!initialize()) {
      return false;
    }

    const int width = framebuffer_width();
    const int height = framebuffer_height();
    if (!valid_backend_dimensions(width, height, actual_backend(), limits_)) {
      report_error("Display framebuffer dimensions exceed configured limits");
      return false;
    }

    if (width != dm_get_width(dmp_) || height != dm_get_height(dmp_)) {
      active_view_->gv_width = width;
      active_view_->gv_height = height;
      dm_set_width(dmp_, width);
      dm_set_height(dmp_, height);
      if (dm_configure_win(dmp_, 0) != BRLCAD_OK) {
        report_error("Unable to resize the libdm framebuffer");
        return false;
      }
      dirty_ = true;
    }

    if (actual_backend() == Backend::OSMesa && !dirty_ && !pixels_.empty()) {
      return true;
    }

    dm_loadmatrix(dmp_, active_view_->gv_model2view, 0);
    if (dm_draw_begin(dmp_) != BRLCAD_OK) {
      report_error("libdm failed to begin a frame");
      return false;
    }
    bool drawn = false;
    try {
      EditPreview preview(ged_, active_view_);
      dm_draw_objs(active_view_, draw_handler_ ? draw_custom : nullptr,
                   draw_handler_ ? this : nullptr);
      if (draw_handler_) {
        dm_hud_begin(dmp_);
        dm_draw_faceplate(active_view_);
        dm_hud_end(dmp_);
      }
      drawn = true;
    } catch (const std::exception& error) {
      report_error(std::string("Faceplate drawing callback failed: ") + error.what());
    } catch (...) {
      report_error("Faceplate drawing callback failed");
    }
    if (!drawn) dm_hud_end(dmp_);
    if (dm_draw_end(dmp_) != BRLCAD_OK) {
      report_error("libdm failed to finish a frame");
      return false;
    }
    if (!drawn) return false;

    if (actual_backend() == Backend::OSMesa && !capture_pixels()) {
      return false;
    }
    dm_set_dirty(dmp_, 0);
    dirty_ = false;
    return true;
  }

  bool capture_pixels()
  {
    const int width = dmp_ ? dm_get_width(dmp_) : 0;
    const int height = dmp_ ? dm_get_height(dmp_) : 0;
    if (!valid_backend_dimensions(width, height, actual_backend(), limits_)) {
      report_error("Display framebuffer dimensions exceed configured limits");
      return false;
    }
    unsigned char* image = nullptr;
    if (!dmp_ || dm_get_display_image(dmp_, &image, 1, 1) != BRLCAD_OK ||
        !image) {
      report_error("Unable to capture the display framebuffer");
      return false;
    }
    const std::size_t bytes = static_cast<std::size_t>(width) * height * 4U;
    pixels_.assign(image, image + bytes);
    bu_free(image, "libfltkcad display image");
    pixel_width_ = width;
    pixel_height_ = height;
    return true;
  }

  std::optional<Image> capture()
  {
    // FLTK controls presentation and may discard the previous back buffer.
    // Capture a fresh frame before it is presented, after applying any resize.
    if (!render() || (gl_canvas_ && !capture_pixels())) {
      return std::nullopt;
    }
    return Image{pixel_width_, pixel_height_, pixels_};
  }

  void set_view(struct bview* borrowed)
  {
    struct bview* replacement = borrowed ? borrowed : local_view_;
    if (replacement == active_view_) {
      return;
    }
    if (active_view_ && active_view_->dmp == dmp_) {
      active_view_->dmp = nullptr;
    }
    active_view_ = replacement;
    if (dmp_) {
      active_view_->dmp = dmp_;
      active_view_->gv_width = framebuffer_width();
      active_view_->gv_height = framebuffer_height();
      dm_set_vp(dmp_, &active_view_->gv_scale);
    }
    invalidate();
  }

  void invalidate()
  {
    dirty_ = true;
    if (dmp_) {
      dm_set_dirty(dmp_, 1);
    }
    canvas_->redraw();
  }

  int handle_input(int event)
  {
    if (!active_view_) {
      return 0;
    }

    const int event_x = Fl::event_x() - canvas_->x();
    const int event_y = Fl::event_y() - canvas_->y();
    const Fl_Cad_View::InputEvent input{
      event, event_x, event_y, Fl::event_button(), Fl::event_buttons(),
      Fl::event_state(), Fl::event_dx(), Fl::event_dy(), Fl::event_key()};
    return process_input(input);
  }

  int process_input(const Fl_Cad_View::InputEvent& input)
  {
    if (!active_view_) return 0;
    if (input.type == FL_PUSH) canvas_->take_focus();
    if (faceplate_gui_) {
      try {
        if (faceplate_input(input)) return 1;
      } catch (const std::exception& error) {
        faceplate_capture_ = 0;
        report_error(std::string("Faceplate input callback failed: ") + error.what());
        return 1;
      } catch (...) {
        faceplate_capture_ = 0;
        report_error("Faceplate input callback failed");
        return 1;
      }
    }
    if (input_handler_) {
      try {
        if (input_handler_(input)) return 1;
      } catch (const std::exception& error) {
        report_error(std::string("Display input callback failed: ") + error.what());
        return 1;
      } catch (...) {
        report_error("Display input callback failed");
        return 1;
      }
    }
    if (!default_bindings_) {
      return 0;
    }

    switch (input.type) {
      case FL_FOCUS:
        return 1;
      case FL_PUSH:
        drag_x_ = input.x;
        drag_y_ = input.y;
        drag_button_ = input.button;
        drag_state_ = input.state & (FL_SHIFT | FL_CTRL | FL_ALT);
        if (drag_state_ & (FL_SHIFT | FL_CTRL)) return 1;
        if (input.button == FL_LEFT_MOUSE) return zoom(kMouseZoomOut);
        if (input.button == FL_RIGHT_MOUSE) return zoom(kMouseZoomIn);
        return 0;
      case FL_RELEASE:
      case FL_UNFOCUS:
        stop_drag();
        return 1;
      case FL_KEYUP:
        if ((input.state & (FL_SHIFT | FL_CTRL | FL_ALT)) != drag_state_) stop_drag();
        return 1;
      case FL_DRAG:
        return drag(input);
      case FL_MOUSEWHEEL:
        return zoom(std::pow(input.wheel_y < 0 ? kWheelZoomIn : kWheelZoomOut,
                             std::abs(input.wheel_y)));
      case FL_KEYDOWN:
        return key(input.key);
      default:
        return 0;
    }
  }

  void stop_drag()
  {
    drag_x_ = drag_y_ = -1;
    drag_button_ = 0;
    drag_state_ = 0;
  }

  int zoom(double factor)
  {
    active_view_->gv_scale = std::max(BV_MINVIEWSCALE, active_view_->gv_scale / factor);
    active_view_->gv_size = 2.0 * active_view_->gv_scale;
    active_view_->gv_isize = 1.0 / active_view_->gv_size;
    bv_update(active_view_);
    changed();
    return 1;
  }

  int drag(const Fl_Cad_View::InputEvent& input)
  {
    if (drag_x_ < 0 || drag_y_ < 0 || drag_button_ < FL_LEFT_MOUSE ||
        drag_button_ > FL_RIGHT_MOUSE ||
        (input.buttons & FL_BUTTON(drag_button_)) == 0) return 0;
    if ((input.state & (FL_SHIFT | FL_CTRL | FL_ALT)) != drag_state_) {
      stop_drag();
      return 1;
    }
    const int dx = input.x - drag_x_;
    const int dy = input.y - drag_y_;
    drag_x_ = input.x;
    drag_y_ = input.y;
    const bool shift = (drag_state_ & FL_SHIFT) != 0;
    const bool control = (drag_state_ & FL_CTRL) != 0;
    const int axis = (drag_state_ & FL_ALT) ? drag_button_ - FL_LEFT_MOUSE : -1;
    if (!shift && !control) return 1;
    if (shift && control) {
      const int delta = axis == Y ? -dy : axis >= 0 ? dx : dx - dy;
      return zoom(std::exp(static_cast<double>(delta) / kScalePixels));
    }
    point_t center;
    MAT_DELTAS_GET_NEG(center, active_view_->gv_center);
    if (axis < 0) {
      if (bv_adjust(active_view_, dx, dy, center, 0, shift ? BV_TRANS : BV_ROT) > 0)
        changed();
      return 1;
    }
    if (shift) {
      vect_t translation = VINIT_ZERO;
      translation[axis] = static_cast<double>(axis == Y ? -dy : dx) *
                             2.0 / std::max(1, active_view_->gv_width);
      vect_t model_translation;
      MAT4X3VEC(model_translation, active_view_->gv_view2model, translation);
      VSUB2(center, center, model_translation);
      MAT_DELTAS_VEC_NEG(active_view_->gv_center, center);
    } else {
      vect_t angles = VINIT_ZERO;
      angles[axis] = static_cast<double>(axis == X ? dy : dx) * kRotationDegreesPerPixel;
      mat_t rotation;
      bn_mat_angles(rotation, angles[X], angles[Y], angles[Z]);
      bn_mat_mul2(rotation, active_view_->gv_rotation);
    }
    bv_update(active_view_);
    changed();
    return 1;
  }

  bool faceplate_input(const Fl_Cad_View::InputEvent& input)
  {
    if (input.type == FL_UNFOCUS) faceplate_capture_ = 0;
    if (input.type != FL_PUSH && input.type != FL_DRAG && input.type != FL_RELEASE)
      return false;
    if (!faceplate_capture_ && (input.type != FL_PUSH ||
        input.button != FL_MIDDLE_MOUSE || (input.state & (FL_SHIFT | FL_CTRL | FL_ALT))))
      return false;
    if (!initialize()) return false;
    const double x = 2.0 * input.x / std::max(1, owner_.w()) - 1.0;
    const double y = 1.0 - 2.0 * input.y / std::max(1, owner_.h());
    if (input.type == FL_PUSH) {
      faceplate_capture_ = dm_faceplate_gui_select(dmp_, x, y, 0);
      faceplate_capture_y_ = y;
      return faceplate_capture_ != 0;
    }
    if (faceplate_capture_ == DM_FACEPLATE_SLIDER_HIT && faceplate_gui_->sliders) {
      const auto* sliders = faceplate_gui_->sliders;
      dm_faceplate_controls_select(sliders, std::clamp(x, sliders->left, sliders->right),
                                   faceplate_capture_y_, 1);
      changed();
    } else if (input.type == FL_RELEASE) {
      dm_faceplate_gui_select(dmp_, x, y, 1);
      changed();
    }
    if (input.type == FL_RELEASE) faceplate_capture_ = 0;
    return true;
  }

  int key(int value)
  {
    vect_t aet;
    bool handled = true;
    switch (value) {
      case '2': VSET(aet, 35.0, -25.0, 0.0); break;
      case '3': VSET(aet, 35.0, 25.0, 0.0); break;
      case '4': VSET(aet, 45.0, 45.0, 0.0); break;
      case 'f':
      case 'F': VSET(aet, 0.0, 0.0, 0.0); break;
      case 't':
      case 'T': VSET(aet, 270.0, 90.0, 0.0); break;
      case 'l':
      case 'L': VSET(aet, 90.0, 0.0, 0.0); break;
      case 'r':
      case 'R': VSET(aet, 270.0, 0.0, 0.0); break;
      default: handled = false; break;
    }
    if (!handled) {
      return 0;
    }
    bv_view_set_aet(active_view_, aet);
    changed();
    return 1;
  }

  void changed()
  {
    invalidate();
    if (changed_handler_) {
      try {
        changed_handler_();
      } catch (const std::exception& error) {
        report_error(std::string("Display change callback failed: ") +
                     error.what());
      } catch (...) {
        report_error("Display change callback failed");
      }
    }
  }

  void report_error(std::string message)
  {
    if (message == last_error_) {
      return;
    }
    last_error_ = std::move(message);
    if (error_handler_) {
      try {
        error_handler_(last_error_);
      } catch (...) {
        // A diagnostic callback cannot be allowed to unwind through FLTK.
      }
    }
  }

  void draw_error(const Fl_Widget& widget)
  {
    if (gl_canvas_ == &widget) {
      const float background = static_cast<float>(kErrorBackground) / 255.0F;
      glClearColor(background, background, background, 1.0F);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      return;
    }
    fl_color(fl_rgb_color(kErrorBackground, kErrorBackground, kErrorBackground));
    fl_rectf(widget.x(), widget.y(), widget.w(), widget.h());
    fl_color(FL_RED);
    fl_draw(last_error_.c_str(), widget.x() + 10, widget.y() + 10,
            std::max(0, widget.w() - 20), std::max(0, widget.h() - 20),
            FL_ALIGN_CENTER | FL_ALIGN_WRAP);
  }

  int framebuffer_width() const
  {
    return gl_canvas_ ? const_cast<GlCanvas*>(gl_canvas_)->framebuffer_width()
                      : canvas_->w();
  }

  int framebuffer_height() const
  {
    return gl_canvas_ ? const_cast<GlCanvas*>(gl_canvas_)->framebuffer_height()
                      : canvas_->h();
  }

  static void draw_custom(struct bview* view, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.draw_handler_) {
      return;
    }
    try {
      self.draw_handler_(view);
    } catch (const std::exception& error) {
      self.report_error(std::string("Display drawing callback failed: ") + error.what());
    } catch (...) {
      self.report_error("Display drawing callback failed");
    }
  }

  static int extgl_make_current(void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.gl_canvas_) {
      return -1;
    }
    self.gl_canvas_->make_host_current();
    return 0;
  }

  static int extgl_get_size(void* data, int* width, int* height)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.gl_canvas_ || !width || !height) {
      return -1;
    }
    *width = self.gl_canvas_->framebuffer_width();
    *height = self.gl_canvas_->framebuffer_height();
    return valid_dimensions(*width, *height, self.limits_) ? 0 : -1;
  }

  static int extgl_capture(void* data, unsigned char** image, int flip, int alpha)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!image || !self.gl_canvas_) {
      return -1;
    }
    const int width = self.gl_canvas_->framebuffer_width();
    const int height = self.gl_canvas_->framebuffer_height();
    if (!valid_dimensions(width, height, self.limits_)) {
      return -1;
    }
    const int channels = alpha ? 4 : 3;
    const std::size_t bytes = static_cast<std::size_t>(width) * height * channels;
    auto* pixels = static_cast<unsigned char*>(bu_malloc(bytes, "extgl capture"));
    for (int attempt = 0;
         attempt < kMaxGlErrorDrain && glGetError() != GL_NO_ERROR;
         ++attempt) {
    }
    GLint read_buffer = GL_BACK;
    GLint pack_alignment = 4;
    glGetIntegerv(GL_READ_BUFFER, &read_buffer);
    glGetIntegerv(GL_PACK_ALIGNMENT, &pack_alignment);
    glReadBuffer(GL_BACK);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, alpha ? GL_RGBA : GL_RGB,
                 GL_UNSIGNED_BYTE, pixels);
    const GLenum error = glGetError();
    glReadBuffer(static_cast<GLenum>(read_buffer));
    glPixelStorei(GL_PACK_ALIGNMENT, pack_alignment);
    if (error != GL_NO_ERROR) {
      bu_free(pixels, "extgl capture");
      return -1;
    }
    if (flip) {
      flip_rows(pixels, width, height, channels);
    }
    *image = pixels;
    return 0;
  }

  static int extgl_draw_string(void* data, const char* text, int size)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.gl_canvas_ || !text) {
      return -1;
    }
    const int scaled_size = std::max(
      1, static_cast<int>(std::lround(size * self.gl_canvas_->pixels_per_unit())));
    gl_font(FL_HELVETICA, scaled_size);
    gl_draw(text);
    return 0;
  }

  static int extgl_string_bbox(void* data, const char* text, int size,
                               double* width, double* height)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.gl_canvas_ || !text || !width || !height) {
      return -1;
    }
    const int scaled_size = std::max(
      1, static_cast<int>(std::lround(size * self.gl_canvas_->pixels_per_unit())));
    gl_font(FL_HELVETICA, scaled_size);
    *width = gl_width(text);
    *height = gl_height();
    return 0;
  }

  Fl_Cad_View& owner_;
  Backend preference_;
  Fl_Widget* canvas_ = nullptr;
  GlCanvas* gl_canvas_ = nullptr;
  struct bview* local_view_ = nullptr;
  struct ged* ged_ = nullptr;
  struct bview* active_view_ = nullptr;
  struct dm* dmp_ = nullptr;
  const struct dm_faceplate_gui* faceplate_gui_ = nullptr;
  int faceplate_capture_ = 0;
  double faceplate_capture_y_ = 0;
  struct dm_extgl_context extgl_context_{};
  Limits limits_;
  Color background_first_{};
  Color background_second_{};
  DrawHandler draw_handler_;
  InputHandler input_handler_;
  ChangedHandler changed_handler_;
  ErrorHandler error_handler_;
  std::string last_error_;
  std::vector<unsigned char> pixels_;
  int pixel_width_ = 0;
  int pixel_height_ = 0;
  bool dirty_ = true;
  bool default_bindings_ = true;
  bool fallback_pending_ = false;
  int drag_x_ = -1;
  int drag_y_ = -1;
  int drag_button_ = 0;
  int drag_state_ = 0;
};

Fl_Cad_View::Fl_Cad_View(int x, int y, int width, int height,
                         Backend preference, const char* label)
  : Fl_Group(x, y, width, height, label),
    impl_(std::make_unique<Impl>(*this, preference))
{
}

void Fl_Cad_View::ged(struct ged* context) { impl_->ged_ = context; invalidate(); }

Fl_Cad_View::~Fl_Cad_View() = default;

bool Fl_Cad_View::initialize() { return impl_->initialize(); }
bool Fl_Cad_View::initialized() const noexcept { return impl_->dmp_ != nullptr; }
Fl_Cad_View::Backend Fl_Cad_View::backend() const noexcept { return impl_->actual_backend(); }
Fl_Cad_View::Backend Fl_Cad_View::preference() const noexcept { return impl_->preference_; }
struct bview* Fl_Cad_View::view() noexcept { return impl_->active_view_; }
const struct bview* Fl_Cad_View::view() const noexcept { return impl_->active_view_; }
void Fl_Cad_View::view(struct bview* value) { impl_->set_view(value); }
struct dm* Fl_Cad_View::display_manager() noexcept { return impl_->dmp_; }
struct fb* Fl_Cad_View::framebuffer() noexcept
{
  return impl_->dmp_ ? dm_get_fb(impl_->dmp_) : nullptr;
}
void Fl_Cad_View::draw_handler(DrawHandler handler) { impl_->draw_handler_ = std::move(handler); }
bool Fl_Cad_View::process_input(const InputEvent& event) { return impl_->process_input(event) != 0; }
void Fl_Cad_View::input_handler(InputHandler handler) { impl_->input_handler_ = std::move(handler); }
void Fl_Cad_View::faceplate_gui(const struct dm_faceplate_gui* gui)
{
  if (impl_->faceplate_gui_ == gui) return;
  impl_->faceplate_gui_ = gui;
  impl_->faceplate_capture_ = 0;
  dm_set_faceplate_gui(impl_->dmp_, gui);
  invalidate();
}
const struct dm_faceplate_gui* Fl_Cad_View::faceplate_gui() const noexcept
{
  return impl_->faceplate_gui_;
}
void Fl_Cad_View::changed_handler(ChangedHandler handler) { impl_->changed_handler_ = std::move(handler); }
void Fl_Cad_View::error_handler(ErrorHandler handler) { impl_->error_handler_ = std::move(handler); }
void Fl_Cad_View::default_bindings(bool enabled) noexcept { impl_->default_bindings_ = enabled; }
bool Fl_Cad_View::default_bindings() const noexcept { return impl_->default_bindings_; }

void Fl_Cad_View::limits(Limits value)
{
  if (value.max_framebuffer_width <= 0 || value.max_framebuffer_height <= 0 ||
      value.max_framebuffer_pixels == 0) {
    impl_->report_error("Display limits must all be greater than zero");
    return;
  }
  impl_->limits_ = value;
  invalidate();
}

Fl_Cad_View::Limits Fl_Cad_View::limits() const noexcept { return impl_->limits_; }
void Fl_Cad_View::background(Color first, Color second)
{
  impl_->background_first_ = first;
  impl_->background_second_ = second;
  if (impl_->dmp_) {
    if (dm_set_bg(impl_->dmp_, first.red, first.green, first.blue, second.red,
                  second.green, second.blue) != BRLCAD_OK) {
      impl_->report_error("Unable to apply the display background color");
    } else {
      impl_->last_error_.clear();
    }
  }
  invalidate();
}
Fl_Cad_View::Color Fl_Cad_View::background_first() const noexcept
{
  return impl_->background_first_;
}
Fl_Cad_View::Color Fl_Cad_View::background_second() const noexcept
{
  return impl_->background_second_;
}
const std::string& Fl_Cad_View::last_error() const noexcept { return impl_->last_error_; }
void Fl_Cad_View::invalidate() { impl_->invalidate(); }
std::optional<Fl_Cad_View::Image> Fl_Cad_View::capture() { return impl_->capture(); }

void Fl_Cad_View::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  impl_->canvas_->resize(x, y, width, height);
  impl_->invalidate();
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
