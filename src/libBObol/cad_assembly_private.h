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

struct BObolCadBatchBuildState {
    SoBRLCadAssembly *assembly = NULL;
    std::vector<Obol::SharedPartUpdate> parts;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>> partIds;
    std::vector<Obol::InstanceUpdate> instances;
    std::vector<Obol::InstanceId> hiddenInstances;
    std::vector<Obol::InstanceId> selectedInstances;
    std::vector<Obol::InstanceId> unpickableInstances;
    int wireCount = 0;
    int shadedCount = 0;
};

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

    void render(SoGLRenderAction *action);
    void getBounds(SoGetBoundingBoxAction *action);
    void pickRay(SoRayPickAction *action);

protected:
    ~SoBRLCadAssembly(void) override;
    SoDetail *createPickDetail(
	const Obol::CadPickDetailRecord &hit) const override;

private:
    friend class SoBRLDatabaseSource;

    struct CompactInstancePresentation {
	std::string payloadKey;
	Obol::PartId activePart = Obol::CadIdBuilder::Root();
	uint8_t channels = 0;
	int activeLevel = -1;
	uint64_t geometryRevision = 0;
	uint64_t appearanceRevision = 0;
	uint64_t placementRevision = 0;
	uint64_t visibilityRevision = 0;
	uint64_t selectionRevision = 0;
    };

    std::map<Obol::InstanceId, InstanceSemantic> semantics;
    SbBool compactPresentationInitialized = FALSE;
    const void *compactPresentationIndex = NULL;
    uint32_t compactPresentationSourceRevision = 0;
    uint32_t compactPresentationInputsRevision = 0;
    uint64_t compactPresentationCadBatchRevision = 0;
    uint64_t compactPresentationPayloadRevision = 0;
    int compactPresentationDrawMode = -1;
    size_t compactWirePresentationCount = 0;
    size_t compactShadedPresentationCount = 0;
    std::unordered_map<Obol::InstanceId, CompactInstancePresentation,
	std::hash<Obol::InstanceId>> compactInstancePresentations;
    std::unordered_map<Obol::PartId, uint8_t, std::hash<Obol::PartId>>
	compactPartChannels;
    std::unordered_set<Obol::PartId, std::hash<Obol::PartId>>
	compactLodParts;
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
