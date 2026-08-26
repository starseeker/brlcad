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

/* Object-lifetime ids are allocated by database_source.cpp.  Keeping the
 * allocation function out of line gives every responsibility-specific source
 * unit one process-wide sequence, while the private Impl constructor remains
 * inline and cannot become part of libBObol's dynamic ABI. */
uint64_t database_source_handle_id(void);

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
    enum class OverviewState {
	Absent,
	Visible,
	RetirementPending,
	Retired
    };

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
	/* Appending/removing occurrences invalidates a scene-wide coverage
	 * proof.  Visibility and in-place geometry/request updates are exact
	 * local mutations and must not make every view rescan the population. */
	SbBool coverageInvalidated = FALSE;
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
    /* reserveCompactOccurrenceCapacity receives the producer-certified leaf
     * count for a streaming epoch.  Index installation also has an exact
     * current count, but only the streaming count excludes the temporary
     * whole-target overview and can therefore prove leaf-frontier coverage. */
    SbBool compactExpectedInstanceCountCertified = FALSE;
    BObolCompactSourceProfile compactSourceProfile;
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
    /* The overview is an early extent cue, not an LoD leaf.  Its explicit
     * lifecycle prevents a retired overview from being mistaken for missing
     * leaf geometry or resurrected by a later presentation overlay. */
    OverviewState compactOverviewState = OverviewState::Absent;
    SbBool compactVisibilityFrontierActive = FALSE;
    SbBool compactVisibilityFrontierDefault = FALSE;
    std::vector<SbString> compactVisibilityFrontier;
    std::vector<SbBool> compactVisibilityFrontierStates;
    std::vector<SbString> compactSelectedPaths;
    std::vector<PresentationOverride> compactPresentationOverrides;
    /* A completed realization transfers its bounded cold-import window to
     * the live source.  Provider-tail materialization claims entries from the
     * stream exactly once, avoiding a timer-dependent lifetime race without
     * making copied occurrence records own full database geometry. */
    std::shared_ptr<BObolCompactOccurrenceStream> compactStagedSourceStream;
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

/* Keep the complete source-private state visible to the responsibility-
 * specific implementation units.  This is still an uninstalled PImpl; the
 * split avoids accessor layers and extra indirection in compact occurrence
 * hot paths. */
struct SoBRLDatabaseSource::Impl :
    BObolCadSourceState,
    BObolCompactOccurrenceRegistryState,
    BObolCadPresentationBridgeState,
    BObolDatabaseSourceSensorState
{
    Impl(void) :
	BObolCadSourceState(
	    database_source_handle_id(), database_source_handle_id())
    {
    }
};

/* Shared implementation helpers used by the realization and inspection
 * units.  They are hidden libBObol symbols, not an installed API. */
SbString source_effective_instance_key(const SoBRLDatabaseSource *source);
SbString source_effective_representation_key(
    const SoBRLDatabaseSource *source);
SbBox3f compact_part_geometry_bounds(
    const std::shared_ptr<const Obol::PartGeometry> &geometry);
SbBox3f database_source_transform_bounds(const SbBox3f &bounds,
    const SbMatrix &matrix);
SbString record_identity_with_revision(const char *identity,
    uint32_t revision);
const SbString &compact_instance_identity(
    const BObolCompactInstanceEntry &entry);
int database_source_float_different(float a, float b);
int database_source_color_equal(const SbColor &a, const SbColor &b);
int database_source_string_equal(const SbString &a, const char *b);
std::string database_source_db_path_without_instance_suffixes(
    const char *path);
const char *database_source_skip_leading_slash(const char *path);
std::string database_source_leaf_component(const SbString &path);
uint64_t compact_next_revision(uint64_t revision);
bool compact_style_equal(const Obol::InstanceStyle &a,
    const Obol::InstanceStyle &b);
Obol::InstanceStyle compact_effective_style(
    const BObolCompactInstanceEntry &entry);
SbBool compact_effective_authored_visibility(
    const BObolCompactInstanceEntry &entry);
SbBool compact_effective_highlight(
    const BObolCompactInstanceEntry &entry);
void compact_sync_shape_summary_state(BObolCompactInstanceEntry &entry);
SbMatrix compact_mesh_asset_matrix(const SoBRLDatabaseSource *source,
    const BObolCompactInstanceEntry &entry);
Obol::InstanceStyle compact_entry_style_from_source(
    const SoBRLDatabaseSource *source,
    const BObolCompactInstanceEntry &entry,
    SbBool selected, SbBool highlighted);
SbBool node_is_source_placement_transform(const SoNode *node);
void realized_vlist_shape_summary(const SoBRLVListShape *shape,
    BObolRealizedShapeSummary &summary);
void realized_mesh_shape_summary(const SoBRLMeshShape *shape,
    BObolRealizedShapeSummary &summary);

#endif /* LIBBOBOL_DATABASE_SOURCE_PRIVATE_H */
