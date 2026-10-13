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

#include "Uri.h"

#include <algorithm>
#include <cctype>
#include <system_error>

namespace fltkcad::detail
{
namespace
{

bool has_scheme(const std::string& uri)
{
  const auto separator = uri.find(':');
  if(separator == std::string::npos || separator == 0 ||
     std::isalpha(static_cast<unsigned char>(uri.front())) == 0)
  {
    return false;
  }
  return std::all_of(uri.begin(), uri.begin() + separator,
                     [](unsigned char ch) {
                       return std::isalnum(ch) != 0 || ch == '+' || ch == '-' ||
                              ch == '.';
                     });
}

bool is_drive_path(const std::string& path)
{
  return path.size() >= 3 &&
         std::isalpha(static_cast<unsigned char>(path[0])) != 0 &&
         path[1] == ':' && (path[2] == '/' || path[2] == '\\');
}

std::string without_fragment_or_query(std::string uri)
{
  const auto suffix = uri.find_first_of("?#");
  if(suffix != std::string::npos)
  {
    uri.erase(suffix);
  }
  return uri;
}

} // namespace

std::optional<std::filesystem::path> local_path(const std::string& uri)
{
  constexpr const char file_scheme[] = "file://";
  if(uri.compare(0, sizeof(file_scheme) - 1, file_scheme) == 0)
  {
    auto path = without_fragment_or_query(uri.substr(sizeof(file_scheme) - 1));
#ifdef _WIN32
    if(path.size() >= 3 && path.front() == '/' && is_drive_path(path.substr(1)))
    {
      path.erase(path.begin());
    }
#endif
    return std::filesystem::path(path);
  }
  if(has_scheme(uri) && !is_drive_path(uri))
  {
    return std::nullopt;
  }
  return std::filesystem::path(without_fragment_or_query(uri));
}

std::string file_uri(const std::filesystem::path& path)
{
  std::error_code error;
  const auto absolute = std::filesystem::absolute(path, error);
  auto uri_path = (error ? path : absolute).lexically_normal().generic_string();
#ifdef _WIN32
  if(is_drive_path(uri_path))
  {
    uri_path.insert(uri_path.begin(), '/');
  }
#endif
  return "file://" + uri_path;
}

std::string resolve_uri(const std::string& reference,
                        const std::string& base_uri)
{
  if(reference.empty())
  {
    return base_uri;
  }
  if(reference.compare(0, 2, "//") == 0)
  {
    const auto separator = base_uri.find(':');
    return separator == std::string::npos
             ? reference
             : base_uri.substr(0, separator + 1) + reference;
  }
  if(has_scheme(reference) || is_drive_path(reference) ||
     std::filesystem::path(reference).is_absolute())
  {
    return reference;
  }
  if(reference.front() == '#')
  {
    return without_fragment_or_query(base_uri) + reference;
  }
  if(reference.front() == '?')
  {
    return without_fragment_or_query(base_uri) + reference;
  }

  const auto base_path = local_path(base_uri);
  if(!base_path)
  {
    const auto slash = base_uri.find_last_of('/');
    return slash == std::string::npos
             ? reference
             : base_uri.substr(0, slash + 1) + reference;
  }

  auto directory = *base_path;
  if(!base_uri.empty() && base_uri.back() != '/')
  {
    directory = directory.parent_path();
  }
  const auto resolved = (directory / reference).lexically_normal();
  return base_uri.compare(0, 7, "file://") == 0 ? file_uri(resolved)
                                                 : resolved.string();
}

} // namespace fltkcad::detail

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
