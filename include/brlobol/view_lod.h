/*                    V I E W _ L O D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/view_lod.h */

#ifndef BRLOBOL_VIEW_LOD_H
#define BRLOBOL_VIEW_LOD_H

#include "brlobol/defines.h"
#include "brlobol/lod_realization.h"

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

/**
 * View-local active LoD bindings for shared Obol geometry.
 *
 * Shared database objects keep their full-detail/source identity in the scene
 * graph.  Each view controller owns one BRLObolViewLodState that records the
 * active display payloads selected for that view.  Coin actions reach the
 * state through SoBRLViewLodElement rather than by mutating shared mesh nodes.
 */
class BRLOBOL_EXPORT BRLObolViewLodState {
public:
    struct BRLOBOL_EXPORT MeshPayload {
	BRLObolLodMeshPayload mesh;
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
	BRLObolLodCounts counts;
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

    BRLObolViewLodState(void);
    ~BRLObolViewLodState(void);

    void clear(void);
    SbBool applyMeshResult(const SoBRLMeshShape *shape,
	    const BRLObolLodResult &result);
    const MeshPayload *findMesh(const SoBRLMeshShape *shape) const;
    const MeshPayload *findMeshForResult(const BRLObolLodResult &result) const;
    size_t bindingCount(void) const;
    size_t payloadCount(void) const;
    size_t estimateDisplayMeshBytes(void) const;
    size_t evictDisplayMeshes(unsigned int *evictedMeshCount = NULL);

private:
    std::unordered_map<std::string, MeshPayloadPtr> meshBindings;
};

class BRLOBOL_EXPORT SoBRLViewLodElement : public SoElement {
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
	    const BRLObolViewLodState *viewState);
    static const BRLObolViewLodState *get(SoState *state);

protected:
    virtual ~SoBRLViewLodElement(void);

private:
    const BRLObolViewLodState *viewState;
};

class BRLOBOL_EXPORT SoBRLViewLodGroup : public SoGroup {
    typedef SoGroup inherited;

    SO_NODE_HEADER(SoBRLViewLodGroup);

public:
    SoBRLViewLodGroup(void);
    static void initClass(void);

    void setViewLodState(BRLObolViewLodState *viewState);
    BRLObolViewLodState *getViewLodState(void) const;

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

    BRLObolViewLodState *viewState;
};

BRLOBOL_EXPORT const BRLObolViewLodState::MeshPayload *
brlobol_view_lod_mesh_for_action(SoAction *action,
	const SoBRLMeshShape *shape);

#endif /* BRLOBOL_VIEW_LOD_H */
