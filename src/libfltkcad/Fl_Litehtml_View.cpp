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

#include <fltkcad/Fl_Litehtml_View.H>

#include "Uri.h"

#include <FL/Fl.H>
#include <FL/Fl_BMP_Image.H>
#include <FL/Fl_GIF_Image.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_PNG_Image.H>
#include <FL/Fl_RGB_Image.H>
#include <FL/Fl_SVG_Image.H>
#include <FL/Fl_Scrollbar.H>
#include <FL/Fl_Window.H>
#include <FL/fl_draw.H>

#include <litehtml.h>
#include <litehtml/render_item.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace
{

constexpr int default_font_size = 16;
constexpr float points_per_inch = 72.0F;
constexpr float screen_dpi = 96.0F;
constexpr int scroll_line_size = 24;
constexpr int scroll_margin = 8;
constexpr unsigned char opaque = 255;
constexpr std::size_t maximum_gradient_pixels = 4U * 1024U * 1024U;
constexpr std::size_t no_match = std::numeric_limits<std::size_t>::max();

Fl_Color all_match_color()
{
  return fl_rgb_color(255, 241, 118);
}

Fl_Color current_match_color()
{
  return fl_rgb_color(255, 174, 66);
}

int rounded(litehtml::pixel_t value)
{
  return static_cast<int>(std::lround(value.value()));
}

int positive_size(litehtml::pixel_t value)
{
  return std::max(0, rounded(value));
}

std::string ascii_lower(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
    return ch < 128 ? static_cast<char>(std::tolower(ch))
                    : static_cast<char>(ch);
  });
  return text;
}

bool is_html_path(const std::filesystem::path& path)
{
  const auto extension = ascii_lower(path.extension().string());
  return extension == ".html" || extension == ".htm" || extension == ".xhtml";
}

std::optional<std::string> read_bounded_file(
  const std::filesystem::path& path, std::size_t byte_limit,
  const std::string& description, std::string& error)
{
  std::error_code status_error;
  const auto status = std::filesystem::status(path, status_error);
  if(status_error)
  {
    error = "Unable to inspect " + description + ": " +
            status_error.message();
    return std::nullopt;
  }
  if(!std::filesystem::is_regular_file(status))
  {
    error = description + " is not a regular file";
    return std::nullopt;
  }

  std::error_code size_error;
  const auto file_size = std::filesystem::file_size(path, size_error);
  if(size_error)
  {
    error = "Unable to determine " + description + " size: " +
            size_error.message();
    return std::nullopt;
  }
  if(file_size > byte_limit ||
     file_size > static_cast<std::uintmax_t>(
                   std::numeric_limits<std::streamsize>::max()))
  {
    error = description + " exceeds byte limit";
    return std::nullopt;
  }

  std::ifstream stream(path, std::ios::binary);
  if(!stream)
  {
    error = "Unable to open " + description;
    return std::nullopt;
  }
  std::string contents(static_cast<std::size_t>(file_size), '\0');
  if(!contents.empty())
  {
    stream.read(contents.data(), static_cast<std::streamsize>(contents.size()));
    if(stream.gcount() != static_cast<std::streamsize>(contents.size()))
    {
      error = "Unable to read " + description;
      return std::nullopt;
    }
  }
  char extra = '\0';
  if(stream.get(extra))
  {
    error = description + " changed while reading or exceeds byte limit";
    return std::nullopt;
  }
  if(!stream.eof())
  {
    error = "Unable to read " + description;
    return std::nullopt;
  }
  return contents;
}

unsigned char composite_over_white(unsigned char channel, unsigned char alpha)
{
  return static_cast<unsigned char>(
    (channel * alpha + opaque * (opaque - alpha)) / opaque);
}

Fl_Color fltk_color(const litehtml::web_color& color)
{
  return fl_rgb_color(composite_over_white(color.red, color.alpha),
                      composite_over_white(color.green, color.alpha),
                      composite_over_white(color.blue, color.alpha));
}

std::string mime_type_for(const std::filesystem::path& path)
{
  const auto extension = ascii_lower(path.extension().string());
  if(extension == ".png")
  {
    return "image/png";
  }
  if(extension == ".jpg" || extension == ".jpeg")
  {
    return "image/jpeg";
  }
  if(extension == ".gif")
  {
    return "image/gif";
  }
  if(extension == ".svg" || extension == ".svgz")
  {
    return "image/svg+xml";
  }
  if(extension == ".bmp")
  {
    return "image/bmp";
  }
  if(extension == ".css")
  {
    return "text/css";
  }
  if(extension == ".html" || extension == ".htm")
  {
    return "text/html";
  }
  return "application/octet-stream";
}

litehtml::web_color interpolate(const litehtml::web_color& first,
                                const litehtml::web_color& second, float amount)
{
  amount = std::clamp(amount, 0.0F, 1.0F);
  const auto channel = [amount](unsigned char a, unsigned char b) {
    return static_cast<unsigned char>(
      std::lround(a + (static_cast<float>(b) - a) * amount));
  };
  return {channel(first.red, second.red), channel(first.green, second.green),
          channel(first.blue, second.blue),
          channel(first.alpha, second.alpha)};
}

litehtml::web_color sample_gradient(
  const std::vector<litehtml::background_layer::color_point>& points,
  float offset)
{
  if(points.empty())
  {
    return litehtml::web_color::transparent;
  }
  if(offset <= points.front().offset)
  {
    return points.front().color;
  }
  for(std::size_t index = 1; index < points.size(); ++index)
  {
    if(offset <= points[index].offset)
    {
      const float span = points[index].offset - points[index - 1].offset;
      const float amount = ZERO(span)
                             ? 0.0F
                             : (offset - points[index - 1].offset) / span;
      return interpolate(points[index - 1].color, points[index].color,
                         amount);
    }
  }
  return points.back().color;
}

} // namespace

class Fl_Litehtml_View::Impl
{
public:
  class Container;

  struct TextRun
  {
    std::string text;
    litehtml::uint_ptr font = 0;
    litehtml::position position;
    std::size_t text_begin = 0;
  };

  struct Match
  {
    std::size_t begin = 0;
    std::size_t end = 0;
  };

  struct Highlight
  {
    std::size_t begin = 0;
    std::size_t end = 0;
    std::size_t match = 0;
  };

  explicit Impl(Fl_Litehtml_View& owner);
  ~Impl();

  bool set_html(std::string source, std::string base_uri);
  bool load_file(const std::string& path);
  void layout();
  bool layout_safely();
  void draw_document();
  int handle(int event);

  std::optional<Fl_Litehtml_View::Resource>
  load_resource(const std::string& uri) const;
  void report_error(const std::string& message) const;
  void anchor_clicked(const std::string& reference);

  void record_text(const char* text, litehtml::uint_ptr font,
                   const litehtml::position& position);
  void collect_text();
  void draw_run_highlights(const char* text, litehtml::uint_ptr font,
                           const litehtml::position& position);
  std::size_t find(const std::string& query,
                   Fl_Litehtml_View::FindOptions options);
  void rebuild_matches(std::size_t preferred_match = 0);
  void rebuild_highlights();
  bool change_match(int direction);
  void show_current_match();
  void clear_find();

  bool scroll_to_fragment(const std::string& fragment);
  void set_scroll(int x, int y);
  litehtml::position viewport() const;

  static void scrollbar_callback(Fl_Widget*, void* data);

  Fl_Litehtml_View& owner_;
  Fl_Scrollbar* vertical_scrollbar_ = nullptr;
  Fl_Scrollbar* horizontal_scrollbar_ = nullptr;
  std::unique_ptr<Container> container_;
  litehtml::document::ptr document_;

  Fl_Litehtml_View::ResourceLoader resource_loader_;
  Fl_Litehtml_View::LinkHandler link_handler_;
  std::optional<std::string> pending_link_;
  Fl_Litehtml_View::ErrorHandler error_handler_;
  Fl_Litehtml_View::Limits limits_;
  mutable std::string last_error_;
  mutable std::size_t error_generation_ = 0;
  mutable std::size_t resource_count_ = 0;
  bool text_collection_truncated_ = false;

  std::string source_;
  std::string base_uri_;
  std::string title_;
  std::string plain_text_;
  std::string query_;
  Fl_Litehtml_View::FindOptions find_options_;

  std::vector<TextRun> text_runs_;
  std::unordered_map<std::string, std::size_t> run_lookup_;
  std::vector<Match> matches_;
  std::vector<std::vector<Highlight>> highlights_;
  std::size_t current_match_ = no_match;

  int viewport_width_ = 0;
  int viewport_height_ = 0;
  int content_width_ = 0;
  int content_height_ = 0;
  int scroll_x_ = 0;
  int scroll_y_ = 0;
};

class Fl_Litehtml_View::Impl::Container final
  : public litehtml::document_container
{
public:
  explicit Container(Impl& impl) : impl_(impl) {}

  struct Font
  {
    Fl_Font face = FL_HELVETICA;
    int size = default_font_size;
    int height = default_font_size;
    int ascent = default_font_size;
    int descent = 0;
    int decoration_line = litehtml::text_decoration_line_none;
    litehtml::web_color decoration_color = litehtml::web_color::current_color;
  };

  litehtml::uint_ptr create_font(const litehtml::font_description& description,
                                 const litehtml::document*,
                                 litehtml::font_metrics* metrics) override
  {
    auto font = std::make_unique<Font>();
    font->size = std::max(1, rounded(description.size));

    const auto family = ascii_lower(description.family);
    const bool monospaced = family.find("mono") != std::string::npos ||
                            family.find("courier") != std::string::npos;
    const bool sans_serif = family.find("sans-serif") != std::string::npos ||
                            family.find("arial") != std::string::npos ||
                            family.find("helvetica") != std::string::npos;
    const bool serif = !monospaced && !sans_serif &&
                       (family.find("serif") != std::string::npos ||
                        family.find("times") != std::string::npos);
    const bool bold = description.weight >= 600;
    const bool italic = description.style == litehtml::font_style_italic;

    if(monospaced)
    {
      font->face = bold && italic ? FL_COURIER_BOLD_ITALIC
                   : bold        ? FL_COURIER_BOLD
                   : italic      ? FL_COURIER_ITALIC
                                 : FL_COURIER;
    }
    else if(serif)
    {
      font->face = bold && italic ? FL_TIMES_BOLD_ITALIC
                   : bold        ? FL_TIMES_BOLD
                   : italic      ? FL_TIMES_ITALIC
                                 : FL_TIMES;
    }
    else
    {
      font->face = bold && italic ? FL_HELVETICA_BOLD_ITALIC
                   : bold        ? FL_HELVETICA_BOLD
                   : italic      ? FL_HELVETICA_ITALIC
                                 : FL_HELVETICA;
    }

    fl_font(font->face, font->size);
    font->height = std::max(1, static_cast<int>(std::ceil(fl_height())));
    font->descent = std::max(0, static_cast<int>(std::ceil(fl_descent())));
    font->ascent = font->height - font->descent;
    font->decoration_line = description.decoration_line;
    font->decoration_color = description.decoration_color;

    metrics->font_size = font->size;
    metrics->height = font->height;
    metrics->ascent = font->ascent;
    metrics->descent = font->descent;
    metrics->x_height = std::max(1, font->height / 2);
    metrics->ch_width = text_width_for(*font, "0");
    metrics->draw_spaces = true;
    metrics->sub_shift = std::max(1, font->height / 5);
    metrics->super_shift = std::max(1, font->height / 3);

    return reinterpret_cast<litehtml::uint_ptr>(font.release());
  }

  void delete_font(litehtml::uint_ptr handle) override
  {
    delete reinterpret_cast<Font*>(handle);
  }

  litehtml::pixel_t text_width(const char* text,
                               litehtml::uint_ptr handle) override
  {
    if(text == nullptr || handle == 0)
    {
      return 0;
    }
    return text_width_for(*reinterpret_cast<Font*>(handle), text);
  }

  float measure(litehtml::uint_ptr handle, std::string_view text)
  {
    if(handle == 0 || text.empty())
    {
      return 0.0F;
    }
    auto* font = reinterpret_cast<Font*>(handle);
    fl_font(font->face, font->size);
    return static_cast<float>(fl_width(text.data(), static_cast<int>(text.size())));
  }

  void draw_text(litehtml::uint_ptr, const char* text,
                 litehtml::uint_ptr handle, litehtml::web_color color,
                 const litehtml::position& position) override
  {
    if(text == nullptr || handle == 0)
    {
      return;
    }
    if(collecting_)
    {
      impl_.record_text(text, handle, position);
      return;
    }

    auto* font = reinterpret_cast<Font*>(handle);
    fl_font(font->face, font->size);
    impl_.draw_run_highlights(text, handle, position);
    fl_color(fltk_color(color));
    const int x = rounded(position.x);
    const int baseline = rounded(position.y) + font->ascent;
    fl_draw(text, x, baseline);

    const int text_width_pixels = positive_size(text_width(text, handle));
    const auto decoration_color = font->decoration_color.is_current_color
                                    ? color
                                    : font->decoration_color;
    fl_color(fltk_color(decoration_color));
    if((font->decoration_line & litehtml::text_decoration_line_underline) != 0)
    {
      fl_line(x, baseline + std::max(1, font->descent / 2),
              x + text_width_pixels, baseline + std::max(1, font->descent / 2));
    }
    if((font->decoration_line & litehtml::text_decoration_line_overline) != 0)
    {
      fl_line(x, rounded(position.y), x + text_width_pixels,
              rounded(position.y));
    }
    if((font->decoration_line &
        litehtml::text_decoration_line_line_through) != 0)
    {
      const int line_y = rounded(position.y) + font->height / 2;
      fl_line(x, line_y, x + text_width_pixels, line_y);
    }
  }

  litehtml::pixel_t pt_to_px(float points) const override
  {
    return points * screen_dpi / points_per_inch;
  }

  litehtml::pixel_t get_default_font_size() const override
  {
    return default_font_size;
  }

  const char* get_default_font_name() const override { return "sans-serif"; }

  void draw_list_marker(litehtml::uint_ptr,
                        const litehtml::list_marker& marker) override
  {
    if(collecting_)
    {
      return;
    }
    const int x = rounded(marker.pos.x);
    const int y = rounded(marker.pos.y);
    const int width = positive_size(marker.pos.width);
    const int height = positive_size(marker.pos.height);
    fl_color(fltk_color(marker.color));
    switch(marker.marker_type)
    {
      case litehtml::list_style_type_circle:
        fl_arc(x, y, width, height, 0, 360);
        return;
      case litehtml::list_style_type_disc:
        fl_pie(x, y, width, height, 0, 360);
        return;
      case litehtml::list_style_type_square:
        fl_rectf(x, y, width, height);
        return;
      case litehtml::list_style_type_none:
        return;
      default:
        break;
    }

    const std::string label = std::to_string(marker.index) + ".";
    auto* font = reinterpret_cast<Font*>(marker.font);
    if(font != nullptr)
    {
      fl_font(font->face, font->size);
    }
    fl_draw(label.c_str(), x, y, width, height, FL_ALIGN_RIGHT | FL_ALIGN_INSIDE);
  }

  void load_image(const char* source, const char* base_url, bool) override
  {
    image(source, base_url);
  }

  void get_image_size(const char* source, const char* base_url,
                      litehtml::size& size) override
  {
    auto* loaded = image(source, base_url);
    if(loaded != nullptr)
    {
      size.width = loaded->w();
      size.height = loaded->h();
    }
  }

  void draw_image(litehtml::uint_ptr, const litehtml::background_layer& layer,
                  const std::string& source,
                  const std::string& base_url) override
  {
    if(collecting_)
    {
      return;
    }
    auto* loaded = image(source.c_str(), base_url.c_str());
    if(loaded == nullptr)
    {
      return;
    }

    int width = positive_size(layer.origin_box.width);
    int height = positive_size(layer.origin_box.height);
    if(width == 0 || height == 0)
    {
      width = loaded->w();
      height = loaded->h();
    }
    if(width <= 0 || height <= 0)
    {
      return;
    }

    Fl_Image* drawable = loaded;
    if(width != loaded->w() || height != loaded->h())
    {
      const std::string key = resolved_uri(source, base_url) + "\n" +
                              std::to_string(width) + "x" +
                              std::to_string(height);
      auto found = scaled_images_.find(key);
      if(found == scaled_images_.end())
      {
        if(images_.size() + scaled_images_.size() >=
           impl_.limits_.max_cached_images)
        {
          impl_.report_error("Scaled image cache limit reached: " + source);
          return;
        }
        const auto scaled_width = static_cast<std::size_t>(width);
        const auto scaled_height = static_cast<std::size_t>(height);
        const std::size_t remaining_pixels =
          decoded_image_pixels_ < impl_.limits_.max_decoded_image_pixels
            ? impl_.limits_.max_decoded_image_pixels - decoded_image_pixels_
            : 0;
        if(scaled_height != 0 && scaled_width > remaining_pixels / scaled_height)
        {
          impl_.report_error("Scaled image cache exceeds pixel limit: " + source);
          return;
        }
        std::unique_ptr<Fl_Image> scaled(loaded->copy(width, height));
        if(!scaled || scaled->fail() != 0)
        {
          impl_.report_error("Unable to scale image: " + source);
          return;
        }
        decoded_image_pixels_ += scaled_width * scaled_height;
        found = scaled_images_.emplace(key, std::move(scaled)).first;
      }
      drawable = found->second.get();
    }
    if(drawable == nullptr)
    {
      return;
    }

    const int clip_x = rounded(layer.clip_box.x);
    const int clip_y = rounded(layer.clip_box.y);
    const int clip_width = positive_size(layer.clip_box.width);
    const int clip_height = positive_size(layer.clip_box.height);
    fl_push_clip(clip_x, clip_y, clip_width, clip_height);

    int start_x = rounded(layer.origin_box.x);
    int start_y = rounded(layer.origin_box.y);
    const bool repeat_x = layer.repeat == litehtml::background_repeat_repeat ||
                          layer.repeat == litehtml::background_repeat_repeat_x;
    const bool repeat_y = layer.repeat == litehtml::background_repeat_repeat ||
                          layer.repeat == litehtml::background_repeat_repeat_y;
    if(repeat_x)
    {
      while(start_x > clip_x)
      {
        start_x -= width;
      }
    }
    if(repeat_y)
    {
      while(start_y > clip_y)
      {
        start_y -= height;
      }
    }

    const int limit_x = repeat_x ? clip_x + clip_width : start_x + 1;
    const int limit_y = repeat_y ? clip_y + clip_height : start_y + 1;
    for(int y = start_y; y < limit_y; y += height)
    {
      for(int x = start_x; x < limit_x; x += width)
      {
        drawable->draw(x, y);
      }
    }
    fl_pop_clip();
  }

  void draw_solid_fill(litehtml::uint_ptr,
                       const litehtml::background_layer& layer,
                       const litehtml::web_color& color) override
  {
    if(collecting_ || color.alpha == 0)
    {
      return;
    }
    fl_push_clip(rounded(layer.clip_box.x), rounded(layer.clip_box.y),
                 positive_size(layer.clip_box.width),
                 positive_size(layer.clip_box.height));
    fl_color(fltk_color(color));
    fl_rectf(rounded(layer.border_box.x), rounded(layer.border_box.y),
             positive_size(layer.border_box.width),
             positive_size(layer.border_box.height));
    fl_pop_clip();
  }

  void draw_linear_gradient(
    litehtml::uint_ptr, const litehtml::background_layer& layer,
    const litehtml::background_layer::linear_gradient& gradient) override
  {
    if(collecting_ || gradient.color_points.empty())
    {
      return;
    }
    const float dx = gradient.end.x - gradient.start.x;
    const float dy = gradient.end.y - gradient.start.y;
    const float denominator = dx * dx + dy * dy;
    draw_gradient_strips(layer, std::abs(dx) >= std::abs(dy),
                         [&](float x, float y) {
                           const float offset = ZERO(denominator)
                                                  ? 0.0F
                                                  : ((x - gradient.start.x) * dx +
                                                     (y - gradient.start.y) * dy) /
                                                      denominator;
                           return sample_gradient(gradient.color_points, offset);
                         });
  }

  void draw_radial_gradient(
    litehtml::uint_ptr, const litehtml::background_layer& layer,
    const litehtml::background_layer::radial_gradient& gradient) override
  {
    if(collecting_ || gradient.color_points.empty())
    {
      return;
    }
    draw_gradient_pixels(layer, [&](float x, float y) {
      const float rx = ZERO(gradient.radius.x) ? 1.0F : gradient.radius.x;
      const float ry = ZERO(gradient.radius.y) ? 1.0F : gradient.radius.y;
      const float nx = (x - gradient.position.x) / rx;
      const float ny = (y - gradient.position.y) / ry;
      return sample_gradient(gradient.color_points, std::sqrt(nx * nx + ny * ny));
    });
  }

  void draw_conic_gradient(
    litehtml::uint_ptr, const litehtml::background_layer& layer,
    const litehtml::background_layer::conic_gradient& gradient) override
  {
    if(collecting_ || gradient.color_points.empty())
    {
      return;
    }
    draw_gradient_pixels(layer, [&](float x, float y) {
      constexpr float pi = 3.14159265358979323846F;
      float angle = std::atan2(x - gradient.position.x,
                               gradient.position.y - y) * 180.0F / pi;
      angle = std::fmod(angle - gradient.angle + 360.0F, 360.0F) / 360.0F;
      return sample_gradient(gradient.color_points, angle);
    });
  }

  void draw_borders(litehtml::uint_ptr, const litehtml::borders& borders,
                    const litehtml::position& position, bool) override
  {
    if(collecting_)
    {
      return;
    }
    draw_horizontal_border(borders.top, rounded(position.x),
                           rounded(position.y), positive_size(position.width));
    draw_horizontal_border(
      borders.bottom, rounded(position.x),
      rounded(position.y + position.height - borders.bottom.width),
      positive_size(position.width));
    draw_vertical_border(borders.left, rounded(position.x),
                         rounded(position.y), positive_size(position.height));
    draw_vertical_border(
      borders.right,
      rounded(position.x + position.width - borders.right.width),
      rounded(position.y), positive_size(position.height));
  }

  void set_caption(const char* caption) override
  {
    impl_.title_ = caption == nullptr ? std::string() : caption;
  }

  void set_base_url(const char* base_url) override
  {
    if(base_url != nullptr)
    {
      document_base_url_ = fltkcad::detail::resolve_uri(base_url,
                                                             impl_.base_uri_);
    }
  }

  void link(const std::shared_ptr<litehtml::document>&,
            const litehtml::element::ptr&) override
  {
  }

  void on_anchor_click(const char* url,
                       const litehtml::element::ptr&) override
  {
    if(url != nullptr)
    {
      impl_.pending_link_ = url;
    }
  }

  void on_mouse_event(const litehtml::element::ptr&,
                      litehtml::mouse_event) override
  {
  }

  void set_cursor(const char* cursor) override
  {
    if(impl_.owner_.window() == nullptr)
    {
      return;
    }
    const std::string name = cursor == nullptr ? std::string() : cursor;
    if(name == "pointer")
    {
      impl_.owner_.window()->cursor(FL_CURSOR_HAND);
    }
    else if(name == "text")
    {
      impl_.owner_.window()->cursor(FL_CURSOR_INSERT);
    }
    else
    {
      impl_.owner_.window()->cursor(FL_CURSOR_DEFAULT);
    }
  }

  void transform_text(std::string& text,
                      litehtml::text_transform transform) override
  {
    if(transform == litehtml::text_transform_none)
    {
      return;
    }
    bool beginning_of_word = true;
    for(char& value : text)
    {
      const auto ch = static_cast<unsigned char>(value);
      if(ch >= 128)
      {
        beginning_of_word = false;
        continue;
      }
      if(transform == litehtml::text_transform_uppercase ||
         (transform == litehtml::text_transform_capitalize && beginning_of_word))
      {
        value = static_cast<char>(std::toupper(ch));
      }
      else if(transform == litehtml::text_transform_lowercase)
      {
        value = static_cast<char>(std::tolower(ch));
      }
      beginning_of_word = std::isspace(ch) != 0;
    }
  }

  void import_css(std::string& text, const std::string& url,
                  std::string& base_url) override
  {
    const auto resolved = resolved_uri(url, base_url);
    const std::size_t previous_error = impl_.error_generation_;
    const auto resource = impl_.load_resource(resolved);
    if(!resource)
    {
      if(impl_.error_generation_ == previous_error)
      {
        impl_.report_error("Unable to load stylesheet: " + resolved);
      }
      return;
    }
    text.assign(reinterpret_cast<const char*>(resource->data.data()),
                resource->data.size());
    base_url = resolved;
  }

  void set_clip(const litehtml::position& position,
                const litehtml::border_radiuses&) override
  {
    if(!collecting_)
    {
      fl_push_clip(rounded(position.x), rounded(position.y),
                   positive_size(position.width), positive_size(position.height));
    }
  }

  void del_clip() override
  {
    if(!collecting_)
    {
      fl_pop_clip();
    }
  }

  void get_viewport(litehtml::position& viewport) const override
  {
    viewport = impl_.viewport();
  }

  litehtml::element::ptr create_element(
    const char*, const litehtml::string_map&,
    const std::shared_ptr<litehtml::document>&) override
  {
    return nullptr;
  }

  void get_media_features(litehtml::media_features& media) const override
  {
    media.type = litehtml::media_type_screen;
    media.width = impl_.viewport_width_;
    media.height = impl_.viewport_height_;
    media.device_width = Fl::w();
    media.device_height = Fl::h();
    media.color = 8;
    media.color_index = 0;
    media.monochrome = 0;
    media.resolution = screen_dpi;
  }

  void get_language(std::string& language, std::string& culture) const override
  {
    language = "en";
    culture.clear();
  }

  void prepare_document(const std::string& base_uri)
  {
    document_base_url_ = base_uri;
    impl_.resource_count_ = 0;
    images_.clear();
    scaled_images_.clear();
    failed_images_.clear();
    decoded_image_pixels_ = 0;
  }

  void begin_collection() { collecting_ = true; }
  void end_collection() { collecting_ = false; }
  void begin_draw(int x, int y)
  {
    drawing_origin_x_ = x;
    drawing_origin_y_ = y;
  }
  int drawing_origin_x() const { return drawing_origin_x_; }
  int drawing_origin_y() const { return drawing_origin_y_; }
  std::string resolve_reference(const std::string& reference) const
  {
    return resolved_uri(reference, {});
  }

private:
  static litehtml::pixel_t text_width_for(const Font& font, const char* text)
  {
    fl_font(font.face, font.size);
    return static_cast<float>(fl_width(text));
  }

  std::string resolved_uri(const std::string& source,
                           const std::string& base_url) const
  {
    const auto& base = base_url.empty() ? document_base_url_ : base_url;
    return fltkcad::detail::resolve_uri(source, base);
  }

  Fl_Image* image(const char* source, const char* base_url)
  {
    if(source == nullptr)
    {
      return nullptr;
    }
    const std::string resolved =
      resolved_uri(source, base_url == nullptr ? std::string() : base_url);
    const auto existing = images_.find(resolved);
    if(existing != images_.end())
    {
      return existing->second.get();
    }
    if(failed_images_.count(resolved) != 0)
    {
      return nullptr;
    }

    if(images_.size() + scaled_images_.size() >=
       impl_.limits_.max_cached_images)
    {
      failed_images_.insert(resolved);
      impl_.report_error("Image cache limit reached while loading: " + resolved);
      return nullptr;
    }

    const std::size_t previous_error = impl_.error_generation_;
    const auto resource = impl_.load_resource(resolved);
    if(!resource || resource->data.empty() ||
       resource->data.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
      failed_images_.insert(resolved);
      if(impl_.error_generation_ == previous_error)
      {
        impl_.report_error("Unable to load image: " + resolved);
      }
      return nullptr;
    }

    const auto mime_type = ascii_lower(resource->mime_type);
    const auto* data = resource->data.data();
    const int length = static_cast<int>(resource->data.size());
    std::unique_ptr<Fl_Image> decoded;
    if(mime_type == "image/png")
    {
      decoded = std::make_unique<Fl_PNG_Image>(resolved.c_str(), data, length);
    }
    else if(mime_type == "image/jpeg" || mime_type == "image/jpg")
    {
      decoded = std::make_unique<Fl_JPEG_Image>(resolved.c_str(), data, length);
    }
    else if(mime_type == "image/gif")
    {
      decoded = std::make_unique<Fl_GIF_Image>(resolved.c_str(), data,
                                               resource->data.size());
    }
    else if(mime_type == "image/svg+xml")
    {
      decoded = std::make_unique<Fl_SVG_Image>(resolved.c_str(), data,
                                               resource->data.size());
    }
    else if(mime_type == "image/bmp")
    {
      decoded = std::make_unique<Fl_BMP_Image>(resolved.c_str(), data, length);
    }

    if(!decoded || decoded->fail() != 0)
    {
      failed_images_.insert(resolved);
      impl_.report_error("Unable to decode image: " + resolved);
      return nullptr;
    }
    const auto decoded_width = static_cast<std::size_t>(std::max(0, decoded->w()));
    const auto decoded_height = static_cast<std::size_t>(std::max(0, decoded->h()));
    const std::size_t remaining_pixels =
      decoded_image_pixels_ < impl_.limits_.max_decoded_image_pixels
        ? impl_.limits_.max_decoded_image_pixels - decoded_image_pixels_
        : 0;
    if(decoded_height != 0 && decoded_width > remaining_pixels / decoded_height)
    {
      failed_images_.insert(resolved);
      impl_.report_error("Decoded image exceeds pixel limit: " + resolved);
      return nullptr;
    }
    decoded_image_pixels_ += decoded_width * decoded_height;
    return images_.emplace(resolved, std::move(decoded)).first->second.get();
  }

  template <typename Sampler>
  void draw_gradient_strips(const litehtml::background_layer& layer,
                            bool vertical_strips, Sampler sampler)
  {
    const int x = rounded(layer.border_box.x);
    const int y = rounded(layer.border_box.y);
    const int width = positive_size(layer.border_box.width);
    const int height = positive_size(layer.border_box.height);
    if(width == 0 || height == 0)
    {
      return;
    }
    fl_push_clip(rounded(layer.clip_box.x), rounded(layer.clip_box.y),
                 positive_size(layer.clip_box.width),
                 positive_size(layer.clip_box.height));
    const int count = vertical_strips ? width : height;
    for(int index = 0; index < count; ++index)
    {
      const float sample_x = vertical_strips ? static_cast<float>(x + index) + 0.5F
                                             : static_cast<float>(x) + width / 2.0F;
      const float sample_y = vertical_strips ? static_cast<float>(y) + height / 2.0F
                                             : static_cast<float>(y + index) + 0.5F;
      fl_color(fltk_color(sampler(sample_x, sample_y)));
      if(vertical_strips)
      {
        fl_rectf(x + index, y, 1, height);
      }
      else
      {
        fl_rectf(x, y + index, width, 1);
      }
    }
    fl_pop_clip();
  }

  template <typename Sampler>
  void draw_gradient_pixels(const litehtml::background_layer& layer,
                            Sampler sampler)
  {
    const int border_x = rounded(layer.border_box.x);
    const int border_y = rounded(layer.border_box.y);
    const int clip_x = rounded(layer.clip_box.x);
    const int clip_y = rounded(layer.clip_box.y);
    const int x = std::max(border_x, clip_x);
    const int y = std::max(border_y, clip_y);
    const int right = std::min(border_x + positive_size(layer.border_box.width),
                               clip_x + positive_size(layer.clip_box.width));
    const int bottom = std::min(border_y + positive_size(layer.border_box.height),
                                clip_y + positive_size(layer.clip_box.height));
    const int width = std::max(0, right - x);
    const int height = std::max(0, bottom - y);
    if(width == 0 || height == 0)
    {
      return;
    }
    if(static_cast<std::size_t>(width) >
       maximum_gradient_pixels / static_cast<std::size_t>(height))
    {
      draw_gradient_strips(layer, false, sampler);
      return;
    }
    const std::size_t pixel_count = static_cast<std::size_t>(width) * height;

    constexpr std::size_t channel_count = 3;
    std::vector<unsigned char> pixels(pixel_count * channel_count);
    std::size_t offset = 0;
    for(int row = 0; row < height; ++row)
    {
      for(int column = 0; column < width; ++column)
      {
        const auto color = sampler(static_cast<float>(x + column) + 0.5F,
                                   static_cast<float>(y + row) + 0.5F);
        pixels[offset++] = composite_over_white(color.red, color.alpha);
        pixels[offset++] = composite_over_white(color.green, color.alpha);
        pixels[offset++] = composite_over_white(color.blue, color.alpha);
      }
    }

    Fl_RGB_Image image(pixels.data(), width, height,
                       static_cast<int>(channel_count));
    image.draw(x, y);
  }

  static int line_style(litehtml::border_style style)
  {
    switch(style)
    {
      case litehtml::border_style_dotted:
        return FL_DOT;
      case litehtml::border_style_dashed:
        return FL_DASH;
      default:
        return FL_SOLID;
    }
  }

  static bool visible(const litehtml::border& border)
  {
    return border.width.value() > 0.0F &&
           border.style != litehtml::border_style_none &&
           border.style != litehtml::border_style_hidden &&
           border.color.alpha != 0;
  }

  static void draw_horizontal_border(const litehtml::border& border, int x,
                                     int y, int length)
  {
    if(!visible(border))
    {
      return;
    }
    const int width = std::max(1, rounded(border.width));
    fl_color(fltk_color(border.color));
    if(border.style == litehtml::border_style_dashed ||
       border.style == litehtml::border_style_dotted)
    {
      fl_line_style(line_style(border.style), width);
      fl_line(x, y + width / 2, x + length, y + width / 2);
      fl_line_style(FL_SOLID);
    }
    else
    {
      fl_rectf(x, y, length, width);
    }
  }

  static void draw_vertical_border(const litehtml::border& border, int x,
                                   int y, int length)
  {
    if(!visible(border))
    {
      return;
    }
    const int width = std::max(1, rounded(border.width));
    fl_color(fltk_color(border.color));
    if(border.style == litehtml::border_style_dashed ||
       border.style == litehtml::border_style_dotted)
    {
      fl_line_style(line_style(border.style), width);
      fl_line(x + width / 2, y, x + width / 2, y + length);
      fl_line_style(FL_SOLID);
    }
    else
    {
      fl_rectf(x, y, width, length);
    }
  }

  Impl& impl_;
  bool collecting_ = false;
  int drawing_origin_x_ = 0;
  int drawing_origin_y_ = 0;
  std::string document_base_url_;
  std::unordered_map<std::string, std::unique_ptr<Fl_Image>> images_;
  std::unordered_map<std::string, std::unique_ptr<Fl_Image>> scaled_images_;
  std::unordered_set<std::string> failed_images_;
  std::size_t decoded_image_pixels_ = 0;
};

Fl_Litehtml_View::Impl::Impl(Fl_Litehtml_View& owner) : owner_(owner)
{
  container_ = std::make_unique<Container>(*this);

  owner_.begin();
  vertical_scrollbar_ = new Fl_Scrollbar(0, 0, 0, 0);
  vertical_scrollbar_->type(FL_VERTICAL);
  vertical_scrollbar_->linesize(scroll_line_size);
  vertical_scrollbar_->callback(scrollbar_callback, this);
  horizontal_scrollbar_ = new Fl_Scrollbar(0, 0, 0, 0);
  horizontal_scrollbar_->type(FL_HORIZONTAL);
  horizontal_scrollbar_->linesize(scroll_line_size);
  horizontal_scrollbar_->callback(scrollbar_callback, this);
  owner_.end();
  vertical_scrollbar_->hide();
  horizontal_scrollbar_->hide();
}

Fl_Litehtml_View::Impl::~Impl()
{
  document_.reset();
}

void Fl_Litehtml_View::Impl::report_error(const std::string& message) const
{
  ++error_generation_;
  if(last_error_ == message)
  {
    return;
  }
  last_error_ = message;
  if(error_handler_)
  {
    try
    {
      error_handler_(message);
    }
    catch(...)
    {
      // Error notification must remain safe when invoked by FLTK callbacks.
    }
  }
}

std::optional<Fl_Litehtml_View::Resource>
Fl_Litehtml_View::Impl::load_resource(const std::string& uri) const
{
  if(resource_count_ >= limits_.max_resource_count)
  {
    report_error("Resource count limit reached while loading: " + uri);
    return std::nullopt;
  }
  ++resource_count_;

  if(resource_loader_)
  {
    std::optional<Resource> resource;
    try
    {
      resource = resource_loader_(uri);
    }
    catch(const std::exception& exception)
    {
      report_error("Resource loader failed for " + uri + ": " +
                   exception.what());
      return std::nullopt;
    }
    catch(...)
    {
      report_error("Resource loader failed for " + uri);
      return std::nullopt;
    }
    if(resource && resource->data.size() > limits_.max_resource_bytes)
    {
      report_error("Resource exceeds byte limit: " + uri);
      return std::nullopt;
    }
    return resource;
  }

  const auto path = fltkcad::detail::local_path(uri);
  if(!path)
  {
    return std::nullopt;
  }
  std::string read_error;
  const auto contents = read_bounded_file(*path, limits_.max_resource_bytes,
                                          "resource", read_error);
  if(!contents)
  {
    report_error(read_error + ": " + uri);
    return std::nullopt;
  }

  Resource resource;
  resource.data.assign(contents->begin(), contents->end());
  resource.mime_type = mime_type_for(*path);
  return resource;
}

bool Fl_Litehtml_View::Impl::set_html(std::string source,
                                      std::string base_uri)
{
  last_error_.clear();
  if(source.size() > limits_.max_document_bytes)
  {
    report_error("HTML document exceeds byte limit");
    owner_.redraw();
    return false;
  }
  if(source.find('\0') != std::string::npos)
  {
    report_error("HTML document contains binary data");
    owner_.redraw();
    return false;
  }
  if(base_uri.empty())
  {
    std::error_code error;
    const auto working_directory = std::filesystem::current_path(error);
    if(error)
    {
      report_error("Unable to determine the working directory: " +
                   error.message());
      return false;
    }
    base_uri = fltkcad::detail::file_uri(working_directory) + "/";
  }

  document_.reset();
  source_ = std::move(source);
  base_uri_ = std::move(base_uri);
  title_.clear();
  plain_text_.clear();
  text_runs_.clear();
  matches_.clear();
  highlights_.clear();
  current_match_ = no_match;
  query_.clear();
  container_->prepare_document(base_uri_);

  try
  {
    document_ = litehtml::document::createFromString(source_, container_.get());
  }
  catch(const std::exception& exception)
  {
    report_error(std::string("Unable to parse HTML: ") + exception.what());
    layout_safely();
    owner_.redraw();
    return false;
  }
  catch(...)
  {
    report_error("Unable to parse HTML");
    layout_safely();
    owner_.redraw();
    return false;
  }

  if(!document_)
  {
    report_error("Unable to parse HTML");
    layout_safely();
    owner_.redraw();
    return false;
  }

  scroll_x_ = 0;
  scroll_y_ = 0;
  if(!layout_safely())
  {
    return false;
  }
  owner_.redraw();
  return true;
}

bool Fl_Litehtml_View::Impl::load_file(const std::string& path)
{
  last_error_.clear();
  const std::filesystem::path file_path(path);
  if(!is_html_path(file_path))
  {
    report_error("Unsupported HTML file type: " + path);
    return false;
  }

  std::string read_error;
  auto source = read_bounded_file(file_path, limits_.max_document_bytes,
                                  "HTML file", read_error);
  if(!source)
  {
    report_error(read_error + ": " + path);
    return false;
  }
  return set_html(std::move(*source),
                  fltkcad::detail::file_uri(std::filesystem::path(path)));
}

void Fl_Litehtml_View::Impl::layout()
{
  const int scrollbar_size = std::max(1, Fl::scrollbar_size());
  bool show_vertical = false;
  bool show_horizontal = false;

  if(!document_)
  {
    viewport_width_ = owner_.w();
    viewport_height_ = owner_.h();
    content_width_ = 0;
    content_height_ = 0;
    vertical_scrollbar_->hide();
    horizontal_scrollbar_->hide();
    return;
  }

  for(int pass = 0; pass < 3; ++pass)
  {
    viewport_width_ = std::max(1, owner_.w() - (show_vertical ? scrollbar_size : 0));
    viewport_height_ = std::max(1, owner_.h() - (show_horizontal ? scrollbar_size : 0));
    document_->render(viewport_width_);
    const auto bounded_extent = [this](double value, const char* dimension) {
      if(!std::isfinite(value) || value < 0.0)
      {
        report_error(std::string("Invalid document ") + dimension);
        return 0;
      }
      if(value > limits_.max_layout_extent)
      {
        report_error(std::string("Document ") + dimension +
                     " exceeds the configured layout limit");
        return limits_.max_layout_extent;
      }
      return static_cast<int>(value);
    };
    content_width_ = bounded_extent(std::ceil(document_->width().value()),
                                    "width");
    content_height_ = bounded_extent(std::ceil(document_->height().value()),
                                     "height");

    const bool needs_horizontal = content_width_ > viewport_width_;
    const bool needs_vertical = content_height_ > viewport_height_;
    if(needs_horizontal == show_horizontal && needs_vertical == show_vertical)
    {
      break;
    }
    show_horizontal = needs_horizontal;
    show_vertical = needs_vertical;
  }

  vertical_scrollbar_->resize(owner_.x() + viewport_width_, owner_.y(),
                              scrollbar_size, viewport_height_);
  horizontal_scrollbar_->resize(owner_.x(), owner_.y() + viewport_height_,
                                viewport_width_, scrollbar_size);
  show_vertical ? vertical_scrollbar_->show() : vertical_scrollbar_->hide();
  show_horizontal ? horizontal_scrollbar_->show() : horizontal_scrollbar_->hide();

  set_scroll(scroll_x_, scroll_y_);
  collect_text();
  rebuild_matches(current_match_ == no_match ? 0 : current_match_);
}

bool Fl_Litehtml_View::Impl::layout_safely()
{
  try
  {
    layout();
    return true;
  }
  catch(const std::exception& exception)
  {
    document_.reset();
    report_error(std::string("Unable to lay out HTML: ") + exception.what());
  }
  catch(...)
  {
    document_.reset();
    report_error("Unable to lay out HTML");
  }
  layout();
  owner_.redraw();
  return false;
}

void Fl_Litehtml_View::Impl::set_scroll(int x, int y)
{
  scroll_x_ = std::clamp(x, 0, std::max(0, content_width_ - viewport_width_));
  scroll_y_ = std::clamp(y, 0, std::max(0, content_height_ - viewport_height_));
  vertical_scrollbar_->value(scroll_y_, viewport_height_, 0,
                             std::max(content_height_, viewport_height_));
  horizontal_scrollbar_->value(scroll_x_, viewport_width_, 0,
                               std::max(content_width_, viewport_width_));
  owner_.redraw();
}

litehtml::position Fl_Litehtml_View::Impl::viewport() const
{
  return {scroll_x_, scroll_y_, viewport_width_, viewport_height_};
}

void Fl_Litehtml_View::Impl::draw_document()
{
  if(!document_)
  {
    return;
  }
  const int draw_x = owner_.x() - scroll_x_;
  const int draw_y = owner_.y() - scroll_y_;
  const litehtml::position clip(owner_.x(), owner_.y(), viewport_width_,
                                viewport_height_);
  fl_push_clip(owner_.x(), owner_.y(), viewport_width_, viewport_height_);
  container_->begin_draw(draw_x, draw_y);
  try
  {
    document_->draw(0, draw_x, draw_y, &clip);
  }
  catch(...)
  {
    fl_pop_clip();
    throw;
  }
  fl_pop_clip();
}

void Fl_Litehtml_View::Impl::scrollbar_callback(Fl_Widget*, void* data)
{
  auto* self = static_cast<Impl*>(data);
  self->set_scroll(self->horizontal_scrollbar_->value(),
                   self->vertical_scrollbar_->value());
}

int Fl_Litehtml_View::Impl::handle(int event)
{
  const bool inside = Fl::event_x() >= owner_.x() &&
                      Fl::event_x() < owner_.x() + viewport_width_ &&
                      Fl::event_y() >= owner_.y() &&
                      Fl::event_y() < owner_.y() + viewport_height_;
  const int document_x = Fl::event_x() - owner_.x() + scroll_x_;
  const int document_y = Fl::event_y() - owner_.y() + scroll_y_;
  const int client_x = Fl::event_x() - owner_.x();
  const int client_y = Fl::event_y() - owner_.y();
  const auto redraw = [this](const litehtml::position&) { owner_.redraw(); };

  switch(event)
  {
    case FL_ENTER:
    case FL_MOVE:
      if(document_ && inside)
      {
        document_->on_mouse_over(document_x, document_y, client_x, client_y,
                                 redraw);
        return 1;
      }
      break;
    case FL_LEAVE:
      if(document_)
      {
        document_->on_mouse_leave(redraw);
        if(owner_.window())
        {
          owner_.window()->cursor(FL_CURSOR_DEFAULT);
        }
        return 1;
      }
      break;
    case FL_PUSH:
      if(document_ && inside && Fl::event_button() == FL_LEFT_MOUSE)
      {
        owner_.take_focus();
        document_->on_lbutton_down(document_x, document_y, client_x, client_y,
                                   redraw);
        return 1;
      }
      break;
    case FL_RELEASE:
      if(document_ && Fl::event_button() == FL_LEFT_MOUSE)
      {
        pending_link_.reset();
        document_->on_lbutton_up(document_x, document_y, client_x, client_y,
                                 redraw);
        // Navigation can replace the document.  Let litehtml finish using its
        // elements and fonts before invoking application callbacks.
        auto reference = std::exchange(pending_link_, std::nullopt);
        if(reference)
        {
          anchor_clicked(*reference);
        }
        return 1;
      }
      break;
    case FL_MOUSEWHEEL:
      set_scroll(scroll_x_ + Fl::event_dx() * scroll_line_size,
                 scroll_y_ + Fl::event_dy() * scroll_line_size);
      return 1;
    case FL_FOCUS:
    case FL_UNFOCUS:
      return 1;
    case FL_KEYDOWN:
      switch(Fl::event_key())
      {
        case FL_Up:
          set_scroll(scroll_x_, scroll_y_ - scroll_line_size);
          return 1;
        case FL_Down:
          set_scroll(scroll_x_, scroll_y_ + scroll_line_size);
          return 1;
        case FL_Left:
          set_scroll(scroll_x_ - scroll_line_size, scroll_y_);
          return 1;
        case FL_Right:
          set_scroll(scroll_x_ + scroll_line_size, scroll_y_);
          return 1;
        case FL_Page_Up:
          set_scroll(scroll_x_, scroll_y_ - viewport_height_);
          return 1;
        case FL_Page_Down:
          set_scroll(scroll_x_, scroll_y_ + viewport_height_);
          return 1;
        case FL_Home:
          set_scroll(0, 0);
          return 1;
        case FL_End:
          set_scroll(scroll_x_, content_height_);
          return 1;
        default:
          break;
      }
      break;
    default:
      break;
  }
  return 0;
}

void Fl_Litehtml_View::Impl::anchor_clicked(const std::string& reference)
{
  if(!reference.empty() && reference.front() == '#')
  {
    scroll_to_fragment(reference.substr(1));
  }
  const std::string resolved =
    container_->resolve_reference(reference);
  if(link_handler_)
  {
    try
    {
      link_handler_(resolved);
    }
    catch(const std::exception& exception)
    {
      report_error(std::string("Link handler failed: ") + exception.what());
    }
    catch(...)
    {
      report_error("Link handler failed");
    }
  }
}

std::string run_key(std::string_view text, litehtml::uint_ptr font, float x,
                    float y)
{
  std::ostringstream key;
  key << font << ':' << std::lround(x * 64.0F) << ':'
      << std::lround(y * 64.0F) << ':' << text;
  return key.str();
}

void Fl_Litehtml_View::Impl::record_text(
  const char* text, litehtml::uint_ptr font,
  const litehtml::position& position)
{
  if(text == nullptr || *text == '\0')
  {
    return;
  }
  const std::size_t text_size = std::strlen(text);
  if(text_runs_.size() >= limits_.max_text_runs ||
     text_size > limits_.max_extracted_text_bytes - plain_text_.size())
  {
    if(!text_collection_truncated_)
    {
      text_collection_truncated_ = true;
      report_error("Extracted document text exceeded the configured limit");
    }
    return;
  }
  TextRun run;
  run.text = text;
  run.font = font;
  run.position = position;
  run.text_begin = plain_text_.size();
  plain_text_ += run.text;
  const std::size_t index = text_runs_.size();
  run_lookup_.emplace(run_key(run.text, font, position.x.value(),
                              position.y.value()),
                      index);
  text_runs_.push_back(std::move(run));
}

void Fl_Litehtml_View::Impl::collect_text()
{
  plain_text_.clear();
  text_runs_.clear();
  run_lookup_.clear();
  text_collection_truncated_ = false;
  if(!document_)
  {
    return;
  }
  container_->begin_collection();
  document_->draw(0, 0, 0, nullptr);
  container_->end_collection();
  highlights_.resize(text_runs_.size());
}

void Fl_Litehtml_View::Impl::draw_run_highlights(
  const char* text, litehtml::uint_ptr font,
  const litehtml::position& position)
{
  const float document_x = position.x.value() - container_->drawing_origin_x();
  const float document_y = position.y.value() - container_->drawing_origin_y();
  auto found = run_lookup_.find(run_key(text, font, document_x, document_y));
  std::size_t run_index = found == run_lookup_.end() ? no_match : found->second;
  if(run_index == no_match)
  {
    constexpr float position_tolerance = 1.0F;
    for(std::size_t index = 0; index < text_runs_.size(); ++index)
    {
      const auto& run = text_runs_[index];
      if(run.font == font && run.text == text &&
         std::abs(run.position.x.value() - document_x) <= position_tolerance &&
         std::abs(run.position.y.value() - document_y) <= position_tolerance)
      {
        run_index = index;
        break;
      }
    }
  }
  if(run_index == no_match || run_index >= highlights_.size())
  {
    return;
  }

  const auto& highlights = highlights_[run_index];
  for(const auto& highlight : highlights)
  {
    const std::string_view whole(text);
    const float before = container_->measure(font, whole.substr(0, highlight.begin));
    const float width = container_->measure(
      font, whole.substr(highlight.begin, highlight.end - highlight.begin));
    const bool current = highlight.match == current_match_;
    fl_color(current ? current_match_color() : all_match_color());
    fl_rectf(rounded(position.x) + static_cast<int>(std::floor(before)),
             rounded(position.y), std::max(1, static_cast<int>(std::ceil(width))),
             positive_size(position.height));
  }
}

std::size_t Fl_Litehtml_View::Impl::find(
  const std::string& query, Fl_Litehtml_View::FindOptions options)
{
  if(query.size() > limits_.max_find_query_bytes)
  {
    query_.clear();
    matches_.clear();
    current_match_ = no_match;
    rebuild_highlights();
    report_error("Find query exceeds the configured byte limit");
    owner_.redraw();
    return 0;
  }
  query_ = query;
  find_options_ = options;
  rebuild_matches();
  owner_.redraw();
  return matches_.size();
}

void Fl_Litehtml_View::Impl::rebuild_matches(std::size_t preferred_match)
{
  matches_.clear();
  current_match_ = no_match;
  if(query_.empty() || plain_text_.empty())
  {
    rebuild_highlights();
    return;
  }

  const std::string haystack =
    find_options_.case_sensitive ? plain_text_ : ascii_lower(plain_text_);
  const std::string needle =
    find_options_.case_sensitive ? query_ : ascii_lower(query_);
  for(std::size_t offset = 0;
      (offset = haystack.find(needle, offset)) != std::string::npos;)
  {
    if(matches_.size() >= limits_.max_find_matches)
    {
      report_error("Find results exceeded the configured match limit");
      break;
    }
    matches_.push_back({offset, offset + needle.size()});
    offset += std::max<std::size_t>(1, needle.size());
  }

  if(!matches_.empty())
  {
    current_match_ = std::min(preferred_match, matches_.size() - 1);
  }
  rebuild_highlights();
  show_current_match();
}

void Fl_Litehtml_View::Impl::rebuild_highlights()
{
  highlights_.assign(text_runs_.size(), {});
  std::size_t first_candidate = 0;
  for(std::size_t match_index = 0; match_index < matches_.size(); ++match_index)
  {
    const auto& match = matches_[match_index];
    while(first_candidate < text_runs_.size() &&
          text_runs_[first_candidate].text_begin +
              text_runs_[first_candidate].text.size() <=
            match.begin)
    {
      ++first_candidate;
    }
    for(std::size_t run_index = first_candidate;
        run_index < text_runs_.size() &&
        text_runs_[run_index].text_begin < match.end;
        ++run_index)
    {
      const auto& run = text_runs_[run_index];
      const std::size_t run_end = run.text_begin + run.text.size();
      const std::size_t begin = std::max(match.begin, run.text_begin);
      const std::size_t end = std::min(match.end, run_end);
      if(begin < end)
      {
        highlights_[run_index].push_back(
          {begin - run.text_begin, end - run.text_begin, match_index});
      }
    }
  }
}

bool Fl_Litehtml_View::Impl::change_match(int direction)
{
  if(matches_.empty())
  {
    return false;
  }
  if(current_match_ == no_match)
  {
    current_match_ = direction < 0 ? matches_.size() - 1 : 0;
  }
  else if(direction > 0)
  {
    if(current_match_ + 1 >= matches_.size())
    {
      if(!find_options_.wrap)
      {
        return false;
      }
      current_match_ = 0;
    }
    else
    {
      ++current_match_;
    }
  }
  else if(current_match_ == 0)
  {
    if(!find_options_.wrap)
    {
      return false;
    }
    current_match_ = matches_.size() - 1;
  }
  else
  {
    --current_match_;
  }
  show_current_match();
  owner_.redraw();
  return true;
}

void Fl_Litehtml_View::Impl::show_current_match()
{
  if(current_match_ == no_match || current_match_ >= matches_.size())
  {
    return;
  }
  const auto& match = matches_[current_match_];
  for(const auto& run : text_runs_)
  {
    const std::size_t run_end = run.text_begin + run.text.size();
    if(match.begin < run_end && match.end > run.text_begin)
    {
      const int top = rounded(run.position.y) - scroll_margin;
      const int bottom = rounded(run.position.y + run.position.height) +
                         scroll_margin;
      int next_y = scroll_y_;
      if(top < scroll_y_)
      {
        next_y = top;
      }
      else if(bottom > scroll_y_ + viewport_height_)
      {
        next_y = bottom - viewport_height_;
      }
      set_scroll(scroll_x_, next_y);
      return;
    }
  }
}

void Fl_Litehtml_View::Impl::clear_find()
{
  query_.clear();
  matches_.clear();
  current_match_ = no_match;
  rebuild_highlights();
  owner_.redraw();
}

bool Fl_Litehtml_View::Impl::scroll_to_fragment(const std::string& fragment)
{
  if(!document_ || fragment.empty())
  {
    return false;
  }

  std::vector<litehtml::element::ptr> pending{document_->root()};
  std::size_t visited = 0;
  while(!pending.empty() && visited < limits_.max_fragment_nodes)
  {
    auto element = std::move(pending.back());
    pending.pop_back();
    ++visited;
    if(!element)
    {
      continue;
    }
    const char* id = element->get_attr("id");
    if(id != nullptr && fragment == id)
    {
      const auto placement = element->get_placement();
      set_scroll(scroll_x_, rounded(placement.y));
      return true;
    }
    for(const auto& child : element->children())
    {
      if(child)
      {
        pending.push_back(child);
      }
    }
  }
  if(!pending.empty())
  {
    report_error("Fragment search exceeded the configured node limit");
  }
  return false;
}

Fl_Litehtml_View::Fl_Litehtml_View(int x, int y, int width, int height,
                                   const char* label)
  : Fl_Group(x, y, width, height, label), impl_(std::make_unique<Impl>(*this))
{
  box(FL_FLAT_BOX);
  color(FL_WHITE);
}

Fl_Litehtml_View::~Fl_Litehtml_View() = default;

bool Fl_Litehtml_View::html(std::string source, std::string base_uri)
{
  return impl_->set_html(std::move(source), std::move(base_uri));
}

bool Fl_Litehtml_View::load(const std::string& path)
{
  return impl_->load_file(path);
}

void Fl_Litehtml_View::resource_loader(ResourceLoader loader)
{
  impl_->resource_loader_ = std::move(loader);
}

void Fl_Litehtml_View::link_handler(LinkHandler handler)
{
  impl_->link_handler_ = std::move(handler);
}

void Fl_Litehtml_View::error_handler(ErrorHandler handler)
{
  impl_->error_handler_ = std::move(handler);
}

void Fl_Litehtml_View::limits(Limits value)
{
  if(value.max_document_bytes == 0 || value.max_resource_bytes == 0 ||
     value.max_resource_count == 0 || value.max_cached_images == 0 ||
     value.max_decoded_image_pixels == 0 || value.max_text_runs == 0 ||
     value.max_extracted_text_bytes == 0 ||
     value.max_find_query_bytes == 0 || value.max_find_matches == 0 ||
     value.max_fragment_nodes == 0 || value.max_layout_extent <= 0)
  {
    impl_->report_error("HTML limits must all be greater than zero");
    return;
  }
  impl_->limits_ = value;
  impl_->container_->prepare_document(impl_->base_uri_);
  if(impl_->source_.size() > value.max_document_bytes)
  {
    impl_->document_.reset();
    impl_->report_error("Current HTML document exceeds the new byte limit");
  }
  impl_->layout_safely();
}

Fl_Litehtml_View::Limits Fl_Litehtml_View::limits() const noexcept
{
  return impl_->limits_;
}

const std::string& Fl_Litehtml_View::last_error() const noexcept
{
  return impl_->last_error_;
}

const std::string& Fl_Litehtml_View::document_title() const
{
  return impl_->title_;
}

const std::string& Fl_Litehtml_View::plain_text() const
{
  return impl_->plain_text_;
}

const std::string& Fl_Litehtml_View::base_uri() const
{
  return impl_->base_uri_;
}

std::size_t Fl_Litehtml_View::find(const std::string& query,
                                   FindOptions options)
{
  return impl_->find(query, options);
}

std::size_t Fl_Litehtml_View::find(const std::string& query)
{
  return impl_->find(query, FindOptions{});
}

bool Fl_Litehtml_View::find_next()
{
  return impl_->change_match(1);
}

bool Fl_Litehtml_View::find_previous()
{
  return impl_->change_match(-1);
}

void Fl_Litehtml_View::clear_find()
{
  impl_->clear_find();
}

std::size_t Fl_Litehtml_View::match_count() const
{
  return impl_->matches_.size();
}

std::size_t Fl_Litehtml_View::current_match() const
{
  return impl_->current_match_ == no_match ? 0 : impl_->current_match_ + 1;
}

bool Fl_Litehtml_View::scroll_to_fragment(const std::string& fragment)
{
  return impl_->scroll_to_fragment(fragment);
}

void Fl_Litehtml_View::draw()
{
  Fl_Group::draw();
  try
  {
    impl_->draw_document();
  }
  catch(const std::exception& exception)
  {
    impl_->report_error(std::string("Unable to draw HTML: ") + exception.what());
  }
  catch(...)
  {
    impl_->report_error("Unable to draw HTML");
  }
}

int Fl_Litehtml_View::handle(int event)
{
  if(Fl_Group::handle(event) != 0)
  {
    return 1;
  }
  try
  {
    return impl_->handle(event);
  }
  catch(const std::exception& exception)
  {
    impl_->report_error(std::string("Unable to handle HTML event: ") +
                        exception.what());
  }
  catch(...)
  {
    impl_->report_error("Unable to handle HTML event");
  }
  return 0;
}

void Fl_Litehtml_View::resize(int x, int y, int width, int height)
{
  Fl_Group::resize(x, y, width, height);
  impl_->layout_safely();
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
