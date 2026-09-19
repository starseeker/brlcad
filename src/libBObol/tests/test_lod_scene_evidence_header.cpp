/* T E S T _ L O D _ S C E N E _ E V I D E N C E _ H E A D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "../lod_scene_evidence_private.h"

static_assert(sizeof(BObolLodCoveragePolicy) > 0,
    "scene evidence policy must be independently compilable");
