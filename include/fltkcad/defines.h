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

#ifndef FLTKCAD_DEFINES_H
#define FLTKCAD_DEFINES_H

#include "common.h"

#ifndef FLTKCAD_EXPORT
#  if defined(FLTKCAD_DLL_EXPORTS) && defined(FLTKCAD_DLL_IMPORTS)
#    error "Only FLTKCAD_DLL_EXPORTS or FLTKCAD_DLL_IMPORTS can be defined, not both."
#  elif defined(FLTKCAD_DLL_EXPORTS)
#    define FLTKCAD_EXPORT COMPILER_DLLEXPORT
#  elif defined(FLTKCAD_DLL_IMPORTS)
#    define FLTKCAD_EXPORT COMPILER_DLLIMPORT
#  else
#    define FLTKCAD_EXPORT
#  endif
#endif

#endif /* FLTKCAD_DEFINES_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
