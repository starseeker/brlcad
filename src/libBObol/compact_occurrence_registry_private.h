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

#include <deque>
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
	instance(Obol::CadIdBuilder::rootInstance()),
	part(),
	localToSource(SbMatrix::identity()),
	geometryTransform(SbMatrix::identity()),
	placementTransform(SbMatrix::identity()),
	localTransform(SbMatrix::identity()),
	authoredVisible(TRUE),
	presentationVisibleValid(FALSE),
	presentationVisible(TRUE),
	visible(TRUE),
	selectable(TRUE),
	selected(FALSE),
	authoredHighlighted(FALSE),
	presentationHighlightedValid(FALSE),
	presentationHighlighted(FALSE),
	highlighted(FALSE),
	presentationTransparencyValid(FALSE),
	presentationTransparency(0.0f),
	wireGeometry(FALSE),
	pointGeometry(FALSE),
	meshGeometry(FALSE),
	viewDependentCsgGeometry(FALSE),
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
    /* The compact occurrence key is assigned before the entry is indexed and
     * is immutable thereafter.  Presentation hot paths retain this value
     * instead of allocating one SbString copy per occurrence. */
    SbString instanceKey;
    /* Authored fields are the reusable source baseline.  Semantic scene
     * presentation rules are overlays and must never destroy that baseline:
     * an endpoint snapshot may replace or clear the overlay after a source
     * has already streamed thousands of occurrences. */
    SbBool authoredVisible;
    SbBool presentationVisibleValid;
    SbBool presentationVisible;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool authoredHighlighted;
    SbBool presentationHighlightedValid;
    SbBool presentationHighlighted;
    SbBool highlighted;
    SbBool presentationTransparencyValid;
    float presentationTransparency;
    SbBool wireGeometry;
    SbBool pointGeometry;
    SbBool meshGeometry;
    SbBool viewDependentCsgGeometry;
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

/*
 * Pointer identity is a useful fast path while the producer-owned immutable
 * geometry is alive, but a raw pointer alone is not an ownership token.  In
 * particular, content-deduplicated occurrence geometry need not be retained
 * by the part library after its occurrence is upgraded.  A later allocation
 * may reuse that address for unrelated geometry.  Pair the lookup key with a
 * weak lifetime witness so an expired/reused address can never inherit the
 * former part identity.
 */
struct BObolCompactGeometryPartIdentity {
    std::weak_ptr<const Obol::PartGeometry> geometry;
    Obol::PartId part;
};

inline bool
bobol_compact_geometry_identity_matches(
    const BObolCompactGeometryPartIdentity &identity,
    const std::shared_ptr<const Obol::PartGeometry> &candidate)
{
    const std::shared_ptr<const Obol::PartGeometry> retained =
	identity.geometry.lock();
    return retained && candidate && retained.get() == candidate.get();
}

inline bool
bobol_compact_geometry_is_resident_progressive(
    const std::shared_ptr<const Obol::PartGeometry> &geometry)
{
    return geometry && geometry->wire && geometry->wire->isProgressive() &&
	geometry->wire->hasProgressiveErrorBounds();
}

struct BObolCompactInstanceIndex {
    BObolCompactInstanceIndex(void) :
	wireCount(0),
	shadedCount(0),
	viewDependentCsgGeometryCount(0),
	sourceMeshRequestCount(0),
	residentProgressiveGeometryCount(0),
	displayLodTargetCount(0)
    {
	sourceBounds.makeEmpty();
    }

    std::map<std::string, Obol::PartId> partIdByKey;
    std::unordered_map<const Obol::PartGeometry *,
	BObolCompactGeometryPartIdentity>
	partIdByGeometry;
    std::vector<BObolCompactPartReference> parts;
    /* The stream supplies its complete target count before merging.  Reserve
     * exactly that population once so a 150k discovery never doubles a live
     * retained instance array on the GUI thread. */
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    /* Entries are deliberately segmented.  Unlike the compact renderer
     * update array, an entry owns strings, semantic state, and optional
     * geometry metadata.  Reserving an entire large target before its first
     * leaf arrives creates a prohibitive one-shot allocation under pressure;
     * a deque preserves entry addresses and lets progressive discovery retain
     * useful partial coverage. */
    std::deque<BObolCompactInstanceEntry> entries;
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
    size_t viewDependentCsgGeometryCount;
    size_t sourceMeshRequestCount;
    size_t residentProgressiveGeometryCount;
    size_t displayLodTargetCount;
};

#endif /* LIBBOBOL_COMPACT_OCCURRENCE_REGISTRY_PRIVATE_H */
