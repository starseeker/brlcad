/*          L O D _ C O O R D I N A T O R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_COORDINATOR_PRIVATE_H
#define LIBBOBOL_LOD_COORDINATOR_PRIVATE_H

#include "common.h"

#include "lod_admission_policy_private.h"
#include "lod_capacity_policy_private.h"
#include "lod_delivery_policy_private.h"
#include "lod_presentation_policy_private.h"
#include "lod_scene_evidence_private.h"
#include "lod_view_policy_private.h"

/* Private provider parameters for the budgeted terminal-display path.  These
 * are request behavior, not persistent asset identity or public API. */
static constexpr const char *BOBOL_LOD_PREPARED_CAD_ONLY_PARAM =
    "display.prepared_cad_only";

#endif /* LIBBOBOL_LOD_COORDINATOR_PRIVATE_H */
