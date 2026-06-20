/*                       D E F I N E S . H
 * BRL-CAD
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
/** @file brlobol/defines.h */

#ifndef BRLOBOL_DEFINES_H
#define BRLOBOL_DEFINES_H

#include "common.h"

#ifndef BRLOBOL_EXPORT
#  if defined(BRLOBOL_DLL_EXPORTS) && defined(BRLOBOL_DLL_IMPORTS)
#    error "Only BRLOBOL_DLL_EXPORTS or BRLOBOL_DLL_IMPORTS can be defined, not both."
#  elif defined(BRLOBOL_DLL_EXPORTS)
#    define BRLOBOL_EXPORT COMPILER_DLLEXPORT
#  elif defined(BRLOBOL_DLL_IMPORTS)
#    define BRLOBOL_EXPORT COMPILER_DLLIMPORT
#  else
#    define BRLOBOL_EXPORT
#  endif
#endif

#endif /* BRLOBOL_DEFINES_H */
