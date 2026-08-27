/* C O M P A C T _ P R E S E N T A T I O N _ S T A G I N G _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_COMPACT_PRESENTATION_STAGING_PRIVATE_H
#define LIBBOBOL_COMPACT_PRESENTATION_STAGING_PRIVATE_H

#include "cad_assembly_private.h"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/* Copy-on-write bookkeeping for one compact assembly publication. */
class BObolCompactPresentationStaging {
public:
    using CompactPresentation =
	SoBRLCadAssembly::CompactInstancePresentation;
    using CompactLayers =
	std::vector<SoBRLCadAssembly::CompactLayerPresentation>;

    BObolCompactPresentationStaging(
	SoBRLCadAssembly *assembly, bool reset);

    uint8_t partChannel(Obol::PartId part) const;
    void setPartChannel(Obol::PartId part, uint8_t channels);
    void addLodPart(Obol::PartId part);

    bool isSpatiallyCulled(Obol::InstanceId instance) const;
    void setSpatiallyCulled(Obol::InstanceId instance, bool culled);

    const CompactPresentation *findPresentation(
	Obol::InstanceId instance) const;
    CompactPresentation &editPresentation(Obol::InstanceId instance);

    const CompactLayers *findLayers(Obol::InstanceId instance) const;
    void setLayers(Obol::InstanceId instance, CompactLayers layers);

    size_t activePartReferences(Obol::PartId part) const;
    void addPartReference(Obol::PartId part);
    void removePartReference(Obol::PartId part);

    std::vector<Obol::InstanceId> expandInstanceSet(
	const std::vector<Obol::InstanceId> &baseSet) const;
    void appendSpatiallyCulled(
	std::vector<Obol::InstanceId> &instances) const;

    void removePart(Obol::PartId part);
    void commit();

private:
    SoBRLCadAssembly *assembly_;
    bool reset_;
    std::unordered_map<Obol::PartId, uint8_t, std::hash<Obol::PartId>>
	partChannels_;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	removedPartChannels_;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> lodParts_;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	removedLodParts_;
    std::unordered_map<Obol::InstanceId, bool,
	std::hash<Obol::InstanceId>> spatialCull_;
    std::unordered_map<Obol::InstanceId, CompactPresentation,
	std::hash<Obol::InstanceId>> presentations_;
    std::unordered_map<Obol::InstanceId, CompactLayers,
	std::hash<Obol::InstanceId>> layerPresentations_;
    std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	removedLayerPresentations_;
    std::unordered_map<Obol::PartId, ptrdiff_t, std::hash<Obol::PartId>>
	activePartReferenceDeltas_;
};

#endif /* LIBBOBOL_COMPACT_PRESENTATION_STAGING_PRIVATE_H */
