/*      C O M P A C T _ O C C U R R E N C E _ R E G I S T R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_COMPACT_OCCURRENCE_REGISTRY_PRIVATE_H
#define LIBBOBOL_COMPACT_OCCURRENCE_REGISTRY_PRIVATE_H

#include "BObol/BDatabaseSource.h"
#include "cad_assembly_private.h"

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

/*
 * Structure-of-arrays-adjacent compact occurrence state.  This registry is a
 * source-private data model, not a public scene-node hierarchy: occurrence
 * records remain contiguous, geometry is shared immutably, and all lookup
 * tables resolve to stable vector indices.
 */
struct BObolCompactInstanceEntry {
    BObolCompactInstanceEntry(void) :
	instance(Obol::CadIdBuilder::Root()),
	part(Obol::CadIdBuilder::Root()),
	localToSource(SbMatrix::identity()),
	geometryTransform(SbMatrix::identity()),
	placementTransform(SbMatrix::identity()),
	localTransform(SbMatrix::identity()),
	authoredVisible(TRUE),
	visible(TRUE),
	selectable(TRUE),
	selected(FALSE),
	highlighted(FALSE),
	wireGeometry(FALSE),
	pointGeometry(FALSE),
	meshGeometry(FALSE),
	lodBacked(FALSE),
	sourceMeshRequestValid(FALSE),
	geometryRevision(1),
	appearanceRevision(1),
	placementRevision(1),
	visibilityRevision(1),
	selectionRevision(1),
	occurrenceIndex(0),
	booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION)
    {
	sourceBounds.makeEmpty();
    }

    Obol::InstanceId instance;
    Obol::PartId part;
    SbMatrix localToSource;
    SbMatrix geometryTransform;
    /* Tree/object placement without geometry normalization or PCA reuse. */
    SbMatrix placementTransform;
    /* Cached geometryTransform * placementTransform composite. */
    SbMatrix localTransform;
    /* Source-local bounds of the retained part after localTransform. */
    SbBox3f sourceBounds;
    SoBRLCadAssembly::InstanceSemantic semantic;
    /*
     * The compact occurrence key is immutable after the entry is indexed.
     * Keep it lazily materialized for the few legacy constructors that only
     * supplied an InstanceId.  Presentation hot paths can then retain a
     * reference instead of allocating one SbString copy per occurrence.
     */
    mutable SbString instanceKey;
    SbBool authoredVisible;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool highlighted;
    SbBool wireGeometry;
    SbBool pointGeometry;
    SbBool meshGeometry;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    BObolSourceMeshRequest sourceMeshRequest;
    std::shared_ptr<const Obol::PartGeometry> geometry;
    BObolRealizedShapeSummary shapeSummary;
    Obol::InstanceStyle normalStyle;
    Obol::InstanceStyle selectedStyle;
    Obol::InstanceStyle highlightedStyle;
    Obol::InstanceStyle style;
    uint64_t geometryRevision;
    uint64_t appearanceRevision;
    uint64_t placementRevision;
    uint64_t visibilityRevision;
    uint64_t selectionRevision;
    uint32_t occurrenceIndex;
    int booleanOperation;
};

struct BObolCompactPartReference {
    Obol::PartId part;
    std::shared_ptr<const Obol::PartGeometry> geometry;
};

struct BObolCompactInstanceIndex {
    BObolCompactInstanceIndex(void) :
	wireCount(0),
	shadedCount(0),
	sourceMeshRequestCount(0)
    {
	sourceBounds.makeEmpty();
    }

    std::map<std::string, Obol::PartId> partIdByKey;
    std::unordered_map<const Obol::PartGeometry *, Obol::PartId>
	partIdByGeometry;
    std::vector<BObolCompactPartReference> parts;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    std::vector<BObolCompactInstanceEntry> entries;
    std::unordered_map<Obol::InstanceId, size_t,
	std::hash<Obol::InstanceId>> entryIndex;
    std::unordered_map<std::string, size_t> entryIndexByKey;
    std::unordered_map<std::string, size_t> entryIndexByPath;
    /*
     * Semantic subtree operations need ordered prefix lookup while compact
     * occurrences are still arriving in arbitrary batches.  Keep the O(1)
     * exact-path table above for merge/upsert and this ordered companion for
     * selection, visibility, highlighting, and edit-frontier ranges.
     */
    std::map<std::string, size_t> entryIndexByOrderedPath;
    std::unordered_map<std::string, std::vector<size_t>> entryIndicesByLeaf;
    std::unordered_map<std::string, std::vector<size_t>>
	entryIndicesBySourceName;
    std::unordered_map<Obol::PartId, size_t, std::hash<Obol::PartId>>
	partReferenceCounts;
    /*
     * Progressive realization is overwhelmingly append-only.  A running
     * box therefore handles the common case in O(1) without the six
     * red-black-tree nodes previously allocated for every occurrence.
     * Replacing a current extremum marks it dirty; the next query rebuilds
     * once from each entry's cached bounds.
     */
    SbBox3f sourceBounds;
    bool sourceBoundsDirty = false;
    int wireCount;
    int shadedCount;
    size_t sourceMeshRequestCount;
};

#endif /* LIBBOBOL_COMPACT_OCCURRENCE_REGISTRY_PRIVATE_H */
