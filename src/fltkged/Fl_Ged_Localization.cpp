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

#include "Fl_Ged_Localization.H"

#include <clocale>

#if defined(FLTKGED_HAVE_GETTEXT)
#  include <cstdlib>
#  include <libintl.h>

namespace {

constexpr const char* kTranslationDomain = "fltkged";

} // namespace
#endif

void fltkged::initialize_locale()
{
  std::setlocale(LC_ALL, "");
  // GED command syntax always uses a period as the decimal separator.
  std::setlocale(LC_NUMERIC, "C");
#if defined(FLTKGED_HAVE_GETTEXT)
  const char* locale_directory = std::getenv("FLTKGED_LOCALE_PATH");
  if (!locale_directory || !*locale_directory) {
#  if defined(FLTKGED_LOCALE_DIRECTORY)
    locale_directory = FLTKGED_LOCALE_DIRECTORY;
#  endif
  }
  if (locale_directory && *locale_directory) {
    bindtextdomain(kTranslationDomain, locale_directory);
  }
  bind_textdomain_codeset(kTranslationDomain, "UTF-8");
  textdomain(kTranslationDomain);
#endif
}

const char* fltkged::translate(const char* message) noexcept
{
  if (!message) {
    return "";
  }
#if defined(FLTKGED_HAVE_GETTEXT)
  return dgettext(kTranslationDomain, message);
#else
  return message;
#endif
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
