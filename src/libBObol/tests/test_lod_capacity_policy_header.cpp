/* T E S T _ L O D _ C A P A C I T Y _ P O L I C Y _ H E A D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "../lod_capacity_policy_private.h"

static_assert(sizeof(BObolLodCapacityEvidence) > 0,
    "capacity policy must be independently compilable");
