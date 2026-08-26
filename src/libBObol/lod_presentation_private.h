/*        L O D _ P R E S E N T A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_PRESENTATION_PRIVATE_H
#define LIBBOBOL_LOD_PRESENTATION_PRIVATE_H

#include "BObol/BLodRealization.h"

/* These helpers select already prepared immutable renderer objects.  They
 * never build geometry and are therefore safe on the presentation owner
 * thread. */
bool bobol_lod_presentation_geometry_supports_cut(
    const std::shared_ptr<const Obol::PartGeometry> &geometry,
    int drawMode, int activeCut);

bool bobol_lod_select_prepared_layers(
    const std::vector<BObolLodPresentationLayer> &available,
    const std::vector<uint32_t> &requiredChunks,
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int drawMode, int activeCut,
    std::vector<BObolLodPresentationLayer> &selected);

#endif /* LIBBOBOL_LOD_PRESENTATION_PRIVATE_H */
