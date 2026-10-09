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

/* Observable stages of one background mesh-LoD producer.  These values are
 * diagnostics only: they neither select work nor form part of cache/request
 * identity.  Keep NONE at zero so default-initialized public status records
 * unambiguously describe the absence of an executing producer. */
enum BObolLodProducerStage {
    BOBOL_LOD_PRODUCER_STAGE_NONE = 0,
    BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP = 1,
    BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION = 2,
    BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING = 3,
    BOBOL_LOD_PRODUCER_STAGE_BOUNDS_ANALYSIS = 4,
    BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION = 5,
    BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION = 6,
    BOBOL_LOD_PRODUCER_STAGE_SPATIAL_CONSTRUCTION = 7,
    BOBOL_LOD_PRODUCER_STAGE_CACHE_PERSISTENCE = 8,
    /* Appended to preserve the numeric values already exposed by the public
     * diagnostic API.  It executes before hashing for an eligible cold
     * source even though its stable value follows the older stages. */
    BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW = 9,
    /* A resident asset serializes cache loading, construction, and prefix
     * publication.  An executing task may wait here before it can inspect
     * that shared asset. */
    BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION = 10,
    BOBOL_LOD_PRODUCER_STAGE_COUNT = 11
};

#endif /* BOBOL_BDEFINES_H */
