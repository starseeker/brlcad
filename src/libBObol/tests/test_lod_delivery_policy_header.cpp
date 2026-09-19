/* T E S T _ L O D _ D E L I V E R Y _ P O L I C Y _ H E A D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "../lod_delivery_policy_private.h"

static_assert(sizeof(BObolLodPresentationTransaction) > 0,
    "delivery policy must be independently compilable");
