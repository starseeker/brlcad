/*             D A T A B A S E _ S O U R C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DATABASE_SOURCE_PRIVATE_H
#define LIBBOBOL_DATABASE_SOURCE_PRIVATE_H

#include "BObol/BDatabaseSource.h"
#include "compact_occurrence_registry_private.h"

#include <deque>
#include <vector>

/*
 * Private source state is separated by responsibility while retaining direct
 * field storage.  SoBRLDatabaseSource::Impl inherits these non-polymorphic
 * components, preserving field access and layout-like locality without adding
 * allocation or indirection to occurrence hot paths.
 */
struct BObolCadSourceState {
    BObolCadSourceState(uint64_t sourceRoutingId,
	uint64_t compactSourceId) :
	routingId(sourceRoutingId),
	compactHandleSourceId(compactSourceId)
    {
    }

    struct db_i *dbip = NULL;
    struct BObolMeshLod *meshLod = NULL;
    /* Immutable object-lifetime identity used only for owner-thread result
     * routing.  Unlike compactHandleSourceId it must not change when a source
     * is re-realized or adopts a detached compact registry. */
    const uint64_t routingId;
    uint64_t compactHandleSourceId;
    SbBool meshLodBoundsValid = FALSE;
    SbVec3f meshLodBoundsMin = SbVec3f(0.0f, 0.0f, 0.0f);
    SbVec3f meshLodBoundsMax = SbVec3f(0.0f, 0.0f, 0.0f);
};

struct BObolCompactOccurrenceRegistryState {
    struct PresentationOverride {
	enum Property {
	    VISIBILITY = 0,
	    HIGHLIGHT,
	    TRANSPARENCY
	};
	Property property = VISIBILITY;
	SbString path;
	BObolCompactPathMatch match = BOBOL_COMPACT_PATH_EXACT;
	SbBool state = FALSE;
	float transparency = 0.0f;
    };

    struct CadBatchDelta {
	uint64_t revision = 0;
	std::vector<size_t> entryIndices;
    };

    struct DisplayMeshLodDelta {
	uint64_t revision = 0;
	std::vector<size_t> entryIndices;
    };

    struct BObolCompactInstanceIndex *compactIndex = NULL;
    struct BObolCompactInstanceIndex *previousCompactIndex = NULL;
    /*
     * Identity of the complete compact population, independent of allocator
     * addresses and the incremental CAD-batch journal.  Index storage may be
     * recycled after erase/redraw and source/inputs revisions may repeat for
     * a newly constructed source.  The presentation bridge consequently
     * treats (source routing id, population epoch) as its reset token.
     */
    uint64_t compactPopulationEpoch = 1;
    size_t compactExpectedInstanceCount = 0;
    uint64_t cadBatchRevision = 1;
    uint64_t cadBatchDeltaFloorRevision = 1;
    size_t cadBatchDeltaEntryCount = 0;
    std::deque<CadBatchDelta> cadBatchDeltas;
    uint64_t displayMeshLodRevision = 1;
    uint64_t displayMeshLodDeltaFloorRevision = 1;
    size_t displayMeshLodDeltaEntryCount = 0;
    std::deque<DisplayMeshLodDelta> displayMeshLodDeltas;
    /*
     * A streamed compact contract is usable before the owning source's full
     * detached realization is complete.  Record the source epoch at the
     * contract handoff so LoD planning can distinguish that legitimate
     * provisional state from retained requests made stale by a later edit.
     */
    SbBool displayMeshLodContractRevisionValid = FALSE;
    uint32_t displayMeshLodContractSourceRevision = 0;
    uint32_t displayMeshLodContractInputsRevision = 0;
    SbBool compactIndexActive = FALSE;
    SbBool compactOccurrenceRegistry = FALSE;
    SbBool compactVisibilityFrontierActive = FALSE;
    SbBool compactVisibilityFrontierDefault = FALSE;
    std::vector<SbString> compactVisibilityFrontier;
    std::vector<SbBool> compactVisibilityFrontierStates;
    std::vector<SbString> compactSelectedPaths;
    std::vector<PresentationOverride> compactPresentationOverrides;
};

struct BObolCadPresentationBridgeState {
    SoBRLCadAssembly *compiledAssembly = NULL;
    SbBool compiledAssemblyDirty = TRUE;
    SbBool compiledAssemblyActive = FALSE;
    SbUniqueId compiledAssemblyNodeId = 0;
    uint64_t compiledCompactStructureSignature = 0;
    uint64_t compiledCompactStyleSignature = 0;
    uint64_t compiledCompactSemanticSignature = 0;
    uint64_t compiledCompactHiddenSignature = 0;
    uint64_t compiledCompactSelectedSignature = 0;
    uint64_t compiledCompactUnpickableSignature = 0;
};

struct BObolDatabaseSourceSensorState {
    SoFieldSensor *pathSensor = NULL;
    SoFieldSensor *instanceKeySensor = NULL;
    SoFieldSensor *representationKeySensor = NULL;
    SoFieldSensor *representationModeSensor = NULL;
    SoFieldSensor *drawModeSensor = NULL;
    SoFieldSensor *tessellationAbsTolSensor = NULL;
    SoFieldSensor *tessellationRelTolSensor = NULL;
    SoFieldSensor *tessellationNormTolSensor = NULL;
    SoFieldSensor *lodBotThresholdSensor = NULL;
    SoFieldSensor *sourceRevisionSensor = NULL;
    SoFieldSensor *inputsRevisionSensor = NULL;
    SoFieldSensor *viewRevisionSensor = NULL;
};

#endif /* LIBBOBOL_DATABASE_SOURCE_PRIVATE_H */
