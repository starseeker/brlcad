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
    enum NormalStyle {
	NORMAL_AUTHORED = 0,
	NORMAL_FLAT = 1,
	NORMAL_SMOOTH = 2
    };

    struct BOBOL_EXPORT MeshPayload {
	BObolLodMeshPayload mesh;
	BObolLodProgressiveMeshPtr progressiveMesh;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString cacheIdentity;
	SbString cacheKey;
	uint64_t sourceContentHash;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int activeLevel;
	int residentLevel;
	int requestedLevel;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbBool shadedCullBackfaces;
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
	BObolLodProgressiveMeshPtr progressiveMesh;
	BObolLodProxy proxy;
	SbString sourcePath;
	SbString sourceName;
	SbString sourceIdentity;
	SbString sourceInstanceKey;
	SbString sourceBindingKey;
	SbString cacheIdentity;
	SbString cacheKey;
	uint64_t sourceContentHash;
	int resultKind;
	int qualityTier;
	int providerStatus;
	int drawMode;
	int activeLevel;
	int residentLevel;
	int requestedLevel;
	uint64_t viewRevision;
	uint64_t policyRevision;
	BObolLodCounts counts;
	SbBox3f bounds;
	SbBool hasSnappedPoints;
	SbBool hasNormals;
	SbBool shadedCullBackfaces;
	SbString diagnostic;

	CadPayload(void);
	SbBool isValid(void) const;
	size_t estimateBytes(void) const;
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
    /* Hot-path variant for scene planning and telemetry.  Compact
     * occurrence identity, rather than container iteration order, provides
     * determinism to those consumers, so avoid sorting thousands of payload
     * strings on every camera epoch. */
    void findCadPayloadsUnordered(const SoBRLDatabaseSource *source,
	std::vector<const CadPayload *> &payloads) const;
    const CadPayload *findCadForOccurrence(
	const SoBRLDatabaseSource *source,
	const SbString &occurrenceKey) const;
    /* Return a retained progressive asset representative for direct
     * occurrence binding.  Asset residency outlives an individual display
     * occurrence, so off-frustum admission does not force a cache reload. */
    const CadPayload *findCadForAsset(
	const SoBRLDatabaseSource *source,
	const SbString &assetPath) const;
    const CadPayload *findCadForResult(const BObolLodResult &result) const;
    const CadPayload *findCadForResult(const SoBRLDatabaseSource *source,
	const BObolLodResult &result) const;
    /* Change only a view-local PoP cut.  No source/cache work is performed;
     * the call succeeds only when the retained progressive asset already
     * contains the requested prefix. */
    SbBool retargetMeshPayload(const MeshPayload *payload, int activeLevel,
	int requestedLevel, uint64_t viewRevision, uint64_t policyRevision);
    SbBool retargetCadPayload(const CadPayload *payload, int activeLevel,
	int requestedLevel, uint64_t viewRevision, uint64_t policyRevision);
    /* Remove one view-local display binding while retaining its shared asset.
     * The source occurrence's structural fallback becomes visible again.
     * Used by scene-budget/frustum admission when an insignificant occurrence
     * should cost zero triangles rather than its minimum populated PoP cut. */
    SbBool removeCadPayload(const CadPayload *payload);
    SbBool removeMeshPayload(const MeshPayload *payload);
    size_t bindingCount(void) const;
    size_t payloadCount(void) const;
    size_t meshPayloadCount(void) const;
    size_t proxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    size_t cadPayloadCount(void) const;
    size_t cadMeshPayloadCount(void) const;
    size_t cadProxyPayloadCount(int proxyKind = BOBOL_LOD_PROXY_NONE) const;
    /* Estimated triangles submitted by the active view-local cuts.  Shared
     * arrays are counted once per displayed occurrence because render cost
     * follows instances, not storage aliases. */
    size_t activeFaceCount(void) const;
    int maximumActiveProgressiveLevel(void) const;
    /* Apply an O(1)-per-assembly render-only ceiling while the precise
     * occurrence allocator catches up with an interactive view. */
    void setCadPresentationProgressiveLodCeiling(int level) const;
    size_t estimateDisplayMeshBytes(void) const;
    /* Report the richest prefix currently drawn by any view-local occurrence
     * of each retained asset.  This is deliberately activeLevel rather than
     * the unconstrained requestedLevel: the service consumes this snapshot
     * only after submission and worker activity are idle, and uses it to
     * reclaim prefixes above the stable scene-budget cut. */
    void residentMeshDemands(
	std::vector<BObolLodResidentDemand> &demands) const;
    uint64_t cadRevision(void) const;
    /* Compact presentations consume occurrence-local changes without
     * rebuilding every leaf.  fullResync is TRUE when the authoritative
     * source state, rather than the returned keys, must be scanned. */
    void cadOccurrenceChangesSince(const SoBRLDatabaseSource *source,
	uint64_t revision, std::vector<SbString> &occurrenceKeys,
	SbBool &fullResync) const;
    void acknowledgeCadOccurrenceChanges(
	const SoBRLDatabaseSource *source, uint64_t revision) const;
    void noteResidentMeshesChanged(void);
    void setNormalStyle(NormalStyle style, float creaseAngleDegrees = 60.0f);
    NormalStyle getNormalStyle(void) const;
    float getNormalCreaseAngle(void) const;
    /* Drop only mesh/full-detail display data.  Coarse proxy bindings remain
     * resident so memory pressure degrades detail without emptying a scene. */
    size_t evictDisplayMeshPayloads(unsigned int *evictedMeshCount = NULL);
    size_t evictDisplayMeshes(unsigned int *evictedMeshCount = NULL);

    /* Presentation nodes are view-owned, not payload-owned.  A stable node
     * lets an active frame remain drawable while individual LoD occurrences
     * are replaced in place. */
    SoCADAssembly *findCadPresentation(const SoBRLDatabaseSource *source,
	SbString *contentKey = NULL) const;
    void setCadPresentation(const SoBRLDatabaseSource *source,
	SoCADAssembly *assembly, const SbString &contentKey = SbString("")) const;

private:
    struct CadPresentation {
	CadPresentation(void) : assembly(NULL), contentKey("") {}
	SoCADAssembly *assembly;
	SbString contentKey;
    };

    void clearCadPresentations(void) const;
    SbBool applyMeshResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applyProxyResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume);
    SbBool applySourceResultInternal(const SoBRLDatabaseSource *source,
	BObolLodResult &result, SbBool consume);
    std::unordered_map<std::string, MeshPayloadPtr> meshBindings;
    std::unordered_map<std::string, ProxyPayloadPtr> proxyBindings;
    /* One authoritative payload per source/occurrence.  cadBindings retains
     * compatibility aliases for source/path/name lookups, but compact LoD
     * planning and replacement must never scan or deduplicate that alias
     * table. */
    std::unordered_map<std::string,
	std::unordered_map<std::string, CadPayloadPtr> > cadSourceBindings;
    std::unordered_map<std::string,
	std::unordered_map<std::string, CadPayloadPtr> > cadAssetBindings;
    std::unordered_map<std::string, CadPayloadPtr> cadBindings;
    struct CadOccurrenceChange {
	uint64_t revision;
	SbString occurrenceKey;
    };
    mutable std::unordered_map<std::string,
	std::vector<CadOccurrenceChange> > cadOccurrenceChanges;
    uint64_t cadFullResyncRevision;
    uint64_t cadBindingsRevision;
    NormalStyle normalStyle;
    float normalCreaseAngle;
    mutable int cadPresentationProgressiveLodCeiling;
    mutable std::unordered_map<std::string, CadPresentation> cadPresentations;
    void noteCadOccurrenceChanged(const std::string &sourceBindingKey,
	const SbString &occurrenceKey);
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
