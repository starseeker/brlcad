/*                    B V I E W L O D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BViewLod.h */

#ifndef BOBOL_BVIEWLOD_H
#define BOBOL_BVIEWLOD_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/elements/SoElement.h>
#include <Inventor/elements/SoSubElement.h>
#include <Inventor/nodes/SoGroup.h>

#include <memory>
#include <stddef.h>
#include <string>
#include <unordered_map>
#include <vector>

class SoAction;
class SoCallbackAction;
class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoNode;
class SoPickAction;
class SoState;
class SoBRLMeshShape;
class SoBRLDatabaseSource;
class SoCADAssembly;

/**
 * View-local active LoD bindings for shared Obol geometry.
 *
 * Shared database objects keep their full-detail/source identity in the scene
 * graph.  Each view controller owns one BObolViewLodState that records the
 * active display payloads selected for that view.  Coin actions reach the
 * state through SoBRLViewLodElement rather than by mutating shared mesh nodes.
 */
class BOBOL_EXPORT BObolViewLodState
{
public:
    struct BOBOL_EXPORT MeshPayload {
	BObolLodMeshPayload mesh;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString cacheIdentity;
	SbString cacheKey;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int activeLevel;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbString diagnostic;

	MeshPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
	int getTriangleCount(void) const;
	SbBool getTriangleVertexIndices(int triangleIndex,
					int &indexA,
					int &indexB,
					int &indexC) const;
	SbBool getTriangle(int triangleIndex,
			   SbVec3f &a,
			   SbVec3f &b,
			   SbVec3f &c) const;
    };
    typedef std::shared_ptr<MeshPayload> MeshPayloadPtr;

    struct BOBOL_EXPORT ProxyPayload {
	BObolLodProxy proxy;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString cacheIdentity;
	SbString cacheKey;
	int resultKind;
	int qualityTier;
	int providerStatus;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbString diagnostic;

	ProxyPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
    };
    typedef std::shared_ptr<ProxyPayload> ProxyPayloadPtr;

    struct BOBOL_EXPORT CadPayload {
	BObolLodMeshPayload mesh;
	BObolLodProxy proxy;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString sourceInstanceKey;
	SbString sourceBindingKey;
	SbString cacheIdentity;
	SbString cacheKey;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int drawMode;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbString diagnostic;
	mutable SoCADAssembly *assembly;
	mutable SbString assemblyKey;

	CadPayload(void);
	~CadPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
	void clearAssembly(void) const;
    };
    typedef std::shared_ptr<CadPayload> CadPayloadPtr;

    BObolViewLodState(void);
    ~BObolViewLodState(void);

    void clear(void);
    SbBool applyMeshResult(const SoBRLMeshShape *shape,
			   const BObolLodResult &result);
    SbBool applyProxyResult(const SoBRLMeshShape *shape,
			    const BObolLodResult &result);
    SbBool applyDisplayResult(const SoBRLMeshShape *shape,
			      const BObolLodResult &result);
    SbBool applySourceResult(const SoBRLDatabaseSource *source,
			     const BObolLodResult &result);
    SbBool consumeDisplayResult(const SoBRLMeshShape *shape,
	BObolLodResult &result);
    SbBool consumeSourceResult(const SoBRLDatabaseSource *source,
	BObolLodResult &result);
    const MeshPayload *findMesh(const SoBRLMeshShape *shape) const;
    const MeshPayload *findMeshForResult(const BObolLodResult &result) const;
    const ProxyPayload *findProxy(const SoBRLMeshShape *shape) const;
    const ProxyPayload *findProxyForResult(const BObolLodResult &result) const;
    const CadPayload *findCad(const SoBRLDatabaseSource *source) const;
    void findCadPayloads(const SoBRLDatabaseSource *source,
	std::vector<const CadPayload *> &payloads) const;
    const CadPayload *findCadForResult(const BObolLodResult &result) const;
    size_t bindingCount(void) const;
    size_t payloadCount(void) const;
    size_t meshPayloadCount(void) const;
    size_t proxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    size_t cadPayloadCount(void) const;
    size_t cadMeshPayloadCount(void) const;
    size_t cadProxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    size_t estimateDisplayMeshBytes(void) const;
    /* Drop only mesh/full-detail display data.  Coarse proxy bindings remain
     * resident so memory pressure degrades detail without emptying a scene. */
    size_t evictDisplayMeshPayloads(unsigned int *evictedMeshCount = NULL);
    size_t evictDisplayMeshes(unsigned int *evictedMeshCount = NULL);

private:
    SbBool applyMeshResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applyProxyResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applySourceResultInternal(const SoBRLDatabaseSource *source,
	BObolLodResult &result, SbBool consume);
    std::unordered_map<std::string, MeshPayloadPtr> meshBindings;
    std::unordered_map<std::string, ProxyPayloadPtr> proxyBindings;
    std::unordered_map<std::string, CadPayloadPtr> cadBindings;
};

class BOBOL_EXPORT SoBRLViewLodElement : public SoElement
{
    typedef SoElement inherited;

    SO_ELEMENT_HEADER(SoBRLViewLodElement);

public:
    static void initClass(void);

    virtual void init(SoState *state);
    virtual void push(SoState *state);
    virtual SbBool matches(const SoElement *element) const;
    virtual SoElement *copyMatchInfo(void) const;

    static void set(SoState *state,
		    SoNode *node,
		    const BObolViewLodState *viewState);
    static const BObolViewLodState *get(SoState *state);

protected:
    virtual ~SoBRLViewLodElement(void);

private:
    const BObolViewLodState *viewState;
};

class BOBOL_EXPORT SoBRLViewLodGroup : public SoGroup
{
    typedef SoGroup inherited;

    SO_NODE_HEADER(SoBRLViewLodGroup);

public:
    SoBRLViewLodGroup(void);
    static void initClass(void);

    void setViewLodState(BObolViewLodState *viewState);
    BObolViewLodState *getViewLodState(void) const;
    void setSoftwareWireMode(int mode);
    int getSoftwareWireMode(void) const;

    virtual void doAction(SoAction *action) override;
    virtual void GLRender(SoGLRenderAction *action) override;
    virtual void callback(SoCallbackAction *action) override;
    virtual void getBoundingBox(SoGetBoundingBoxAction *action) override;
    virtual void pick(SoPickAction *action) override;

protected:
    virtual ~SoBRLViewLodGroup(void);

private:
    SbBool pushViewState(SoAction *action);
    void popViewState(SoAction *action, SbBool pushed);

    BObolViewLodState *viewState;
    int softwareWireMode;
};

BOBOL_EXPORT const BObolViewLodState::MeshPayload *
bobol_view_lod_mesh_for_action(SoAction *action,
				 const SoBRLMeshShape *shape);

BOBOL_EXPORT const BObolViewLodState::ProxyPayload *
bobol_view_lod_proxy_for_action(SoAction *action,
				  const SoBRLMeshShape *shape);

BOBOL_EXPORT const BObolViewLodState::CadPayload *
bobol_view_lod_cad_for_action(SoAction *action,
				const SoBRLDatabaseSource *source);

BOBOL_EXPORT const BObolViewLodState *
bobol_view_lod_state_for_action(SoAction *action);

#endif /* BOBOL_BVIEWLOD_H */
