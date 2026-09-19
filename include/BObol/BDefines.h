/*                    B D E F I N E S . H
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
/** @file BObol/BDefines.h */

#ifndef BOBOL_BDEFINES_H
#define BOBOL_BDEFINES_H

#include "common.h"

#ifndef BOBOL_EXPORT
#  if defined(BOBOL_DLL_EXPORTS) && defined(BOBOL_DLL_IMPORTS)
#    error "Only BOBOL_DLL_EXPORTS or BOBOL_DLL_IMPORTS can be defined, not both."
#  elif defined(BOBOL_DLL_EXPORTS)
#    define BOBOL_EXPORT COMPILER_DLLEXPORT
#  elif defined(BOBOL_DLL_IMPORTS)
#    define BOBOL_EXPORT COMPILER_DLLIMPORT
#  else
#    define BOBOL_EXPORT
#  endif
#endif

#endif /* BOBOL_BDEFINES_H */
