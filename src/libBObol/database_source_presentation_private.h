/* D A T A B A S E _ S O U R C E _ P R E S E N T A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DATABASE_SOURCE_PRESENTATION_PRIVATE_H
#define LIBBOBOL_DATABASE_SOURCE_PRESENTATION_PRIVATE_H

#include "BObol/BDatabaseSource.h"
#include "BObol/BPickDetail.h"
#include "cad_assembly_private.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec3f.h>

#include <cstddef>
#include <cstdint>
#include <deque>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

class SoAction;

/* Complete owner-thread mutation journal for one compact presentation pass.
 * Publication is derived from these payloads, never from an independently
 * maintained "changed" latch.  Keeping that predicate beside the journal
 * prevents a new sparse mutation category (most notably a cut-only update)
 * from changing BRL-CAD bookkeeping without committing the corresponding
 * Obol transaction. */
struct BObolCadPresentationMutation {
    bool reset = false;
    bool instanceSetsChanged = false;
    bool drawModeChanged = false;
    std::vector<Obol::PartUpdate> baseParts;
    std::vector<Obol::PartUpdate> lodParts;
    std::deque<size_t> changedInstanceIndices;
    std::vector<Obol::InstanceUpdate> layerInstances;
    std::vector<Obol::InstanceId> removedLayerInstances;
    std::vector<Obol::InstanceStyleUpdate> styles;
    std::vector<Obol::InstanceLodUpdate> cuts;
    std::vector<std::pair<Obol::InstanceId,
	SoBRLCadAssembly::InstanceSemantic>> layerSemantics;
    std::vector<Obol::InstanceId> removedLayerSemantics;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> removedParts;

    bool publicationRequired(void) const
    {
	return this->reset || this->instanceSetsChanged ||
	    this->drawModeChanged || !this->baseParts.empty() ||
	    !this->lodParts.empty() ||
	    !this->changedInstanceIndices.empty() ||
	    !this->layerInstances.empty() ||
	    !this->removedLayerInstances.empty() || !this->styles.empty() ||
	    !this->cuts.empty() || !this->layerSemantics.empty() ||
	    !this->removedLayerSemantics.empty() ||
	    !this->removedParts.empty();
    }
};

int source_record_draw_mode(const SoBRLDatabaseSource *source);
std::string cad_instance_key(const SoBRLDatabaseSource *source,
    const char *path, int ordinal);
SbMatrix cad_instance_matrix(const SoBRLDatabaseSource *source,
    const SbMatrix &localMatrix);
Obol::InstanceStyle cad_source_style(const SoBRLDatabaseSource *source);
const char *cad_source_leaf_name(const SoBRLDatabaseSource *source);
Obol::InstanceId cad_source_parent_instance(
    const SoBRLDatabaseSource *source);
uint8_t cad_source_boolean_operation(const SoBRLDatabaseSource *source);
SoBRLCadAssembly::InstanceSemantic cad_source_semantic(
    const SoBRLDatabaseSource *source,
    SoBRLPickDetail::PrimitiveKind primitiveKind);

Obol::CadDrawMode cad_presentation_draw_mode(int sourceDrawMode,
    size_t shadedCount, size_t wireCount);
SoBRLCadAssembly *cad_view_lod_assembly_for_action(SoAction *action,
    const SoBRLDatabaseSource *source);

/* Build the twelve edges of binary-XYZ box corners.  The corner order is the
 * shared aggregate-proxy contract used by Obol validation and PoP metadata. */
int bobol_cad_wire_geometry_from_binary_corners(
    const SbVec3f corners[8], Obol::PartGeometryBuilder &geometry);

#endif /* LIBBOBOL_DATABASE_SOURCE_PRESENTATION_PRIVATE_H */
