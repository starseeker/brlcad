/* T E S T _ L O D _ A D M I S S I O N _ P O L I C Y _ H E A D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "../lod_admission_policy_private.h"

static_assert(sizeof(BObolLodAdmissionPlanner) > 0,
    "admission policy must be independently compilable");
