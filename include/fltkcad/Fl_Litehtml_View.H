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

#ifndef FLTKCAD_FL_LITEHTML_VIEW_H
#define FLTKCAD_FL_LITEHTML_VIEW_H

#include "fltkcad/defines.h"

#include <FL/Fl_Group.H>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class FLTKCAD_EXPORT Fl_Litehtml_View : public Fl_Group
{
public:
  struct Resource
  {
    std::vector<unsigned char> data;
    std::string mime_type;
  };

  struct FindOptions
  {
    bool case_sensitive = false;
    bool wrap = true;
  };

  struct Limits
  {
    std::size_t max_document_bytes = 8U * 1024U * 1024U;
    std::size_t max_resource_bytes = 16U * 1024U * 1024U;
    std::size_t max_resource_count = 256;
    std::size_t max_cached_images = 128;
    std::size_t max_decoded_image_pixels = 32U * 1024U * 1024U;
    std::size_t max_text_runs = 100000;
    std::size_t max_extracted_text_bytes = 16U * 1024U * 1024U;
    std::size_t max_find_query_bytes = 64U * 1024U;
    std::size_t max_find_matches = 100000;
    std::size_t max_fragment_nodes = 250000;
    int max_layout_extent = 1000000;
  };

  using ResourceLoader =
    std::function<std::optional<Resource>(const std::string& uri)>;
  using LinkHandler = std::function<void(const std::string& uri)>;
  using ErrorHandler = std::function<void(const std::string& message)>;

  Fl_Litehtml_View(int x, int y, int width, int height,
                   const char* label = nullptr);
  ~Fl_Litehtml_View() override;

  Fl_Litehtml_View(const Fl_Litehtml_View&) = delete;
  Fl_Litehtml_View& operator=(const Fl_Litehtml_View&) = delete;

  bool html(std::string source, std::string base_uri = {});
  bool load(const std::string& path);

  void resource_loader(ResourceLoader loader);
  void link_handler(LinkHandler handler);
  void error_handler(ErrorHandler handler);

  void limits(Limits value);
  Limits limits() const noexcept;
  const std::string& last_error() const noexcept;

  const std::string& document_title() const;
  const std::string& plain_text() const;
  const std::string& base_uri() const;

  std::size_t find(const std::string& query);
  std::size_t find(const std::string& query, FindOptions options);
  bool find_next();
  bool find_previous();
  void clear_find();
  std::size_t match_count() const;
  std::size_t current_match() const;

  bool scroll_to_fragment(const std::string& fragment);

protected:
  void draw() override;
  int handle(int event) override;
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
