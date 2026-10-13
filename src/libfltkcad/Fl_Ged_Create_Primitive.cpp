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

#include "fltkcad/Fl_Ged_Create_Primitive.H"

#include <FL/Fl_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Input.H>
#include <FL/fl_ask.H>

#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

extern "C" {
#include "bu/vls.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/func.h"
}

namespace {

constexpr int kWindowWidth = 470;
constexpr int kWindowHeight = 205;
constexpr int kMargin = 12;
constexpr int kRowHeight = 28;
constexpr int kLabelWidth = 82;
constexpr int kButtonWidth = 90;
constexpr int kTextLimit = 1024;

bool parse_origin(const char* text, double (&origin)[3])
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  std::string trailing;
  return (input >> origin[0] >> origin[1] >> origin[2]) &&
         !(input >> trailing) && std::isfinite(origin[0]) &&
         std::isfinite(origin[1]) && std::isfinite(origin[2]);
}

bool parse_scale(const char* text, double& scale)
{
  if (!text) {
    return false;
  }
  std::istringstream input(text);
  std::string trailing;
  return (input >> scale) && !(input >> trailing) && std::isfinite(scale) &&
         scale > 0.0;
}

std::string vector_text(const double (&values)[3])
{
  std::ostringstream output;
  output << std::setprecision(12) << values[0] << ' ' << values[1] << ' '
         << values[2];
  return output.str();
}

std::string number_text(double value)
{
  std::ostringstream output;
  output << std::setprecision(12) << value;
  return output.str();
}

} // namespace

class Fl_Ged_Create_Primitive::Impl
{
public:
  explicit Impl(Fl_Ged_Create_Primitive& owner) : owner_(owner)
  {
    owner_.begin();
    name_ = new Fl_Input(kMargin + kLabelWidth, kMargin,
                         owner_.w() - kLabelWidth - 2 * kMargin, kRowHeight,
                         "Name");
    type_ = new Fl_Choice(kMargin + kLabelWidth, kMargin + kRowHeight + 4,
                          owner_.w() - kLabelWidth - 2 * kMargin, kRowHeight,
                          "Type");
    origin_ = new Fl_Input(kMargin + kLabelWidth,
                           kMargin + 2 * (kRowHeight + 4),
                           owner_.w() - kLabelWidth - 2 * kMargin, kRowHeight,
                           "Origin");
    scale_ = new Fl_Input(kMargin + kLabelWidth,
                          kMargin + 3 * (kRowHeight + 4),
                          owner_.w() - kLabelWidth - 2 * kMargin, kRowHeight,
                          "Scale");
    name_->maximum_size(kTextLimit);
    origin_->maximum_size(kTextLimit);
    scale_->maximum_size(kTextLimit);
    origin_->value("0 0 0");
    scale_->value("1");
    create_ = new Fl_Button(owner_.w() - 2 * kButtonWidth - 2 * kMargin,
                            owner_.h() - kRowHeight - kMargin, kButtonWidth,
                            kRowHeight, "Create");
    close_ = new Fl_Button(owner_.w() - kButtonWidth - kMargin,
                           owner_.h() - kRowHeight - kMargin, kButtonWidth,
                           kRowHeight, "Close");
    create_->callback(create_callback, this);
    close_->callback(close_callback, this);
    owner_.end();
    populate_types();
    set_enabled(false);
  }

  void ged(struct ged* value)
  {
    ged_ = value;
    set_enabled(ged_ && ged_->dbip);
    last_error_.clear();
  }

  struct ged* ged() const noexcept { return ged_; }

  bool create()
  {
    if (!ged_ || !ged_->dbip) {
      return fail("Open a database before creating a primitive");
    }
    if (ged_->dbip->dbi_read_only || db_version(ged_->dbip) < 5) {
      return fail("The database is read-only or does not support editing");
    }
    const char* object_name = name_->value();
    const char* primitive_type = type_->text();
    if (!object_name || !*object_name ||
        std::strlen(object_name) > static_cast<std::size_t>(kTextLimit)) {
      return fail("Enter a valid object name");
    }
    if (!primitive_type || !*primitive_type) {
      return fail("Select a primitive type");
    }
    if (db_lookup(ged_->dbip, object_name, LOOKUP_QUIET)) {
      return fail("An object with that name already exists");
    }
    double local_origin[3]{};
    double local_scale = 0.0;
    if (!parse_origin(origin_->value(), local_origin)) {
      return fail("Origin must contain three finite numbers");
    }
    if (!parse_scale(scale_->value(), local_scale)) {
      return fail("Scale must be a positive finite number");
    }
    const double local2base = ged_->dbip->dbi_local2base;
    const double model_origin[3]{local_origin[0] * local2base,
                                 local_origin[1] * local2base,
                                 local_origin[2] * local2base};
    const std::string origin_text = vector_text(model_origin);
    const std::string scale_text = number_text(local_scale * local2base);
    const char* arguments[] = {"make", "-o", origin_text.c_str(), "-s",
                               scale_text.c_str(), object_name, primitive_type};
    bu_vls_trunc(ged_->ged_result_str, 0);
    const int result = ged_exec(ged_, 7, arguments);
    if (result != BRLCAD_OK) {
      const std::string message = bu_vls_strlen(ged_->ged_result_str)
                                    ? bu_vls_cstr(ged_->ged_result_str)
                                    : "BRL-CAD could not create the primitive";
      bu_vls_trunc(ged_->ged_result_str, 0);
      return fail(message);
    }
    bu_vls_trunc(ged_->ged_result_str, 0);
    last_error_.clear();
    const std::string created_name(object_name);
    name_->value("");
    if (created_handler_) {
      created_handler_(created_name);
    }
    return true;
  }

  const std::string& last_error() const noexcept { return last_error_; }
  void created_handler(CreatedHandler handler)
  {
    created_handler_ = std::move(handler);
  }
  void error_handler(ErrorHandler handler)
  {
    error_handler_ = std::move(handler);
  }

private:
  void populate_types()
  {
    struct bu_vls labels = BU_VLS_INIT_ZERO;
    rt_obj_make_labels(&labels, " ");
    std::istringstream input(bu_vls_cstr(&labels));
    std::string label;
    while (input >> label) {
      type_->add(label.c_str());
    }
    bu_vls_free(&labels);
    if (type_->size() > 0) {
      type_->value(0);
    }
  }

  void set_enabled(bool enabled)
  {
    Fl_Widget* widgets[] = {name_, type_, origin_, scale_, create_};
    for (auto* widget : widgets) {
      enabled ? widget->activate() : widget->deactivate();
    }
  }

  bool fail(std::string message)
  {
    last_error_ = std::move(message);
    if (error_handler_) {
      error_handler_(last_error_);
    }
    return false;
  }

  static void create_callback(Fl_Widget*, void* data)
  {
    auto& self = *static_cast<Impl*>(data);
    if (!self.create()) {
      fl_alert("%s", self.last_error_.c_str());
    }
  }
  static void close_callback(Fl_Widget*, void* data)
  {
    static_cast<Impl*>(data)->owner_.hide();
  }

  Fl_Ged_Create_Primitive& owner_;
  Fl_Input* name_ = nullptr;
  Fl_Choice* type_ = nullptr;
  Fl_Input* origin_ = nullptr;
  Fl_Input* scale_ = nullptr;
  Fl_Button* create_ = nullptr;
  Fl_Button* close_ = nullptr;
  struct ged* ged_ = nullptr;
  std::string last_error_;
  CreatedHandler created_handler_;
  ErrorHandler error_handler_;
};

Fl_Ged_Create_Primitive::Fl_Ged_Create_Primitive()
  : Fl_Double_Window(kWindowWidth, kWindowHeight, "Create Primitive"),
    impl_(std::make_unique<Impl>(*this))
{
}
Fl_Ged_Create_Primitive::~Fl_Ged_Create_Primitive() = default;
void Fl_Ged_Create_Primitive::ged(struct ged* value) { impl_->ged(value); }
struct ged* Fl_Ged_Create_Primitive::ged() const noexcept
{
  return impl_->ged();
}
bool Fl_Ged_Create_Primitive::create() { return impl_->create(); }
const std::string& Fl_Ged_Create_Primitive::last_error() const noexcept
{
  return impl_->last_error();
}
void Fl_Ged_Create_Primitive::created_handler(CreatedHandler handler)
{
  impl_->created_handler(std::move(handler));
}
void Fl_Ged_Create_Primitive::error_handler(ErrorHandler handler)
{
  impl_->error_handler(std::move(handler));
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
