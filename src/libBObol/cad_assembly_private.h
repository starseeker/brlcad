/*              C A D _ A S S E M B L Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file cad_assembly_private.h
 *
 * libBObol-private BRL-CAD semantic layer for Obol CAD assemblies.
 */

#ifndef BOBOL_CAD_ASSEMBLY_PRIVATE_H
#define BOBOL_CAD_ASSEMBLY_PRIVATE_H

#include "BObol/BPickDetail.h"

#include <Obol/cad/SoCADAssembly.h>
#include <Obol/cad/CadViewState.h>
#include <Inventor/nodes/SoSeparator.h>

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoRayPickAction;
class SoBRLDatabaseSource;
class BObolViewLodState;
class SoBRLCadAssembly;
class BObolCompactPresentationStaging;

class SoBRLCadAssembly : public SoCADAssembly {
    typedef SoCADAssembly inherited;

    SO_NODE_HEADER(SoBRLCadAssembly);

public:
    struct InstanceSemantic {
	InstanceSemantic(void) :
	    sourceId(0),
	    regionId(0),
	    airCode(0),
	    materialId(0),
	    los(0),
	    materialColorValid(FALSE),
	    primitiveKind(SoBRLPickDetail::UNKNOWN)
	{
	}

	SbString path;
	SbString sourceInstanceKey;
	SbString sourceName;
	SbString sourceType;
	SbString materialShader;
	SbString editIntentId;
	SbString editIntentRole;
	SbColor materialColor;
	uint32_t sourceId;
	int regionId;
	int airCode;
	int materialId;
	int los;
	SbBool materialColorValid;
	SoBRLPickDetail::PrimitiveKind primitiveKind;
    };

    SoBRLCadAssembly(void);
    static void initClass(void);

    void clearSemanticMap(void);
    void setInstanceSemantic(Obol::InstanceId id,
	const InstanceSemantic &semantic);

    /* Exact source occurrence identities for the visible structural repair
     * frontier published by the last complete Obol camera classification. */
    void lastUncollapsedStructuralProxyOccurrenceKeys(
	std::vector<SbString> &occurrenceKeys) const;

    void render(SoGLRenderAction *action);
    void getBounds(SoGetBoundingBoxAction *action);
    void pickRay(SoRayPickAction *action);
    void setPresentationDrawMode(Obol::CadDrawMode mode);
    Obol::CadDrawMode presentationDrawMode(void) const
    {
	return this->presentationDrawModeValue;
    }
    void setPresentationPickMode(Obol::CadPickMode mode);
    Obol::CadPickMode presentationPickMode(void) const
    {
	return this->presentationPickModeValue;
    }

protected:
    ~SoBRLCadAssembly(void) override;
    void GLRender(SoGLRenderAction *action) override;
    void rayPick(SoRayPickAction *action) override;
    SoDetail *createPickDetail(
	const Obol::CadPickDetailRecord &hit) const override;

private:
    friend class SoBRLDatabaseSource;
    friend class BObolCompactPresentationStaging;

    Obol::CadDrawMode presentationDrawModeValue;
    Obol::CadPickMode presentationPickModeValue;

    struct CompactInstancePresentation {
	std::string payloadKey;
	Obol::PartId activePart;
	uint8_t channels = 0;
	int activeCut = -1;
	bool lodStructuralProxy = false;
	uint64_t geometryRevision = 0;
	uint64_t appearanceRevision = 0;
	uint64_t placementRevision = 0;
	uint64_t visibilityRevision = 0;
	uint64_t selectionRevision = 0;
    };

    struct CompactLayerPresentation {
	std::string layerKey;
	Obol::PartId part;
	Obol::InstanceId instance = Obol::CadIdBuilder::rootInstance();
	uint8_t channels = 0;
	int activeCut = -1;
	bool coverage = false;
	uint64_t geometryRevision = 0;
    };

    void reserveCompactPresentationCapacity(size_t expectedOccurrences);

    std::unordered_map<Obol::InstanceId, InstanceSemantic,
	std::hash<Obol::InstanceId>> semantics;
    SbBool compactPresentationInitialized = FALSE;
    const void *compactPresentationIndex = NULL;
    uint64_t compactPresentationSourceRoutingId = 0;
    uint64_t compactPresentationPopulationEpoch = 0;
    uint32_t compactPresentationSourceRevision = 0;
    uint32_t compactPresentationInputsRevision = 0;
    uint64_t compactPresentationCadBatchRevision = 0;
    uint64_t compactPresentationPayloadRevision = 0;
    int compactPresentationDrawMode = -1;
    size_t compactWirePresentationCount = 0;
    size_t compactShadedPresentationCount = 0;
    std::unordered_map<Obol::InstanceId, CompactInstancePresentation,
	std::hash<Obol::InstanceId>> compactInstancePresentations;
    /* Allocated only for the uncommon large-leaf live publication path.  A
     * normal many-leaf scene pays no per-occurrence vector overhead. */
    std::unordered_map<Obol::InstanceId,
	std::vector<CompactLayerPresentation>,
	std::hash<Obol::InstanceId>> compactLayerPresentations;
    /* Exact active presentation references make sparse box->mesh swaps
     * independent of total scene size. */
    std::unordered_map<Obol::PartId, size_t, std::hash<Obol::PartId>>
	compactActivePartReferences;
    std::unordered_map<Obol::PartId, uint8_t, std::hash<Obol::PartId>>
	compactPartChannels;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	compactLodParts;
    /* A spatially paged leaf may have no page intersecting the current view
     * even when its conservative whole-asset bounds touch the frustum.  Keep
     * that view-local cull separate from authored/source visibility so a
     * later camera epoch can restore the occurrence without rediscovery. */
    std::unordered_set<Obol::InstanceId, std::hash<Obol::InstanceId>>
	compactSpatiallyCulledInstances;
};

/** Side-effect-free staging state for one retained assembly publication. */
struct BObolCadBatchBuildState {
    bool valid = true;
    std::vector<Obol::PartUpdate> parts;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> partIds;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<std::pair<Obol::InstanceId,
	SoBRLCadAssembly::InstanceSemantic>> semantics;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    int wireCount = 0;
    int shadedCount = 0;
};

class SoBRLCadRenderBatch : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLCadRenderBatch);

public:
    SoBRLCadRenderBatch(void);
    static void initClass(void);
    void setBatchSourceRoot(SoNode *root);
    void setSoftwareWireMode(int mode);

protected:
    ~SoBRLCadRenderBatch(void) override;
    void GLRender(SoGLRenderAction *action) override;
    void GLRenderBelowPath(SoGLRenderAction *action) override;

private:
    void renderBatch(SoGLRenderAction *action);
    SbBool syncBatch(const BObolViewLodState *viewState);

    SoBRLCadAssembly *assembly;
    SoNode *sourceRoot;
    uint64_t cachedSourceSignature;
    uint64_t cachedStructureSignature;
    uint64_t cachedStyleSignature;
    uint64_t cachedSemanticSignature;
    uint64_t cachedHiddenSignature;
    uint64_t cachedSelectedSignature;
    uint64_t cachedUnpickableSignature;
    SbBool batchValid;
    int softwareWireMode;
    std::unordered_set<const SoBRLDatabaseSource *> batchedSources;
};

int bobol_cad_batch_source_suppressed(const SoBRLDatabaseSource *source);

#endif /* BOBOL_CAD_ASSEMBLY_PRIVATE_H */
