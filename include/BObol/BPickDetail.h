/*                 B P I C K D E T A I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BPickDetail.h */

#ifndef BOBOL_BPICKDETAIL_H
#define BOBOL_BPICKDETAIL_H

#include "BObol/BDefines.h"
#include "BObol/BLodRealization.h"
#include "BObol/BSourceMeshRequest.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>
#include <Inventor/details/SoDetail.h>
#include <Inventor/details/SoSubDetail.h>

#include <stdint.h>
#include <vector>

class BObolLodService;
class SoBRLMeshShape;
struct db_i;
struct resource;
struct rt_i;

class BOBOL_EXPORT SoBRLPickDetail : public SoDetail {
    typedef SoDetail inherited;

    SO_DETAIL_HEADER(SoBRLPickDetail);

public:
    enum PrimitiveKind {
	UNKNOWN = 0,
	LINE_SEGMENT = 1,
	POINT = 2,
	FACE = 3,
	IMPLICIT_SOLID = 4
    };

    SoBRLPickDetail(void);
    SoBRLPickDetail(const SoBRLPickDetail &other);
    SoBRLPickDetail &operator=(const SoBRLPickDetail &other);
    virtual ~SoBRLPickDetail(void);

    static void initClass(void);
    SoDetail *copy(void) const override;
    void clear(void);

    void setPath(const SbString &path);
    const SbString &getPath(void) const;

    void setSourceInstanceKey(const SbString &key);
    const SbString &getSourceInstanceKey(void) const;

    void setSourceName(const SbString &name);
    const SbString &getSourceName(void) const;

    void setSourceType(const SbString &type);
    const SbString &getSourceType(void) const;

    void setSourceId(uint32_t id);
    uint32_t getSourceId(void) const;

    void setRegionId(int id);
    int getRegionId(void) const;

    void setAirCode(int air);
    int getAirCode(void) const;

    void setMaterialId(int material);
    int getMaterialId(void) const;

    void setLos(int los);
    int getLos(void) const;

    void setMaterialColor(SbBool valid, const SbColor &color);
    SbBool hasMaterialColor(void) const;
    const SbColor &getMaterialColor(void) const;

    void setMaterialShader(const SbString &shader);
    const SbString &getMaterialShader(void) const;

    void setPrimitive(PrimitiveKind kind, int index);
    PrimitiveKind getPrimitiveKind(void) const;
    int getPrimitiveIndex(void) const;

    void setEditIntent(const SbString &id, const SbString &role);
    const SbString &getEditIntentId(void) const;
    const SbString &getEditIntentRole(void) const;

    void setFaceVertexIndices(int indexA, int indexB, int indexC);
    int getFaceVertexIndex(int vertexSlot) const;
    int getFaceVertexIndexA(void) const;
    int getFaceVertexIndexB(void) const;
    int getFaceVertexIndexC(void) const;

    void setNearestFaceEdge(int edgeSlot, int indexA, int indexB);
    int getNearestFaceEdgeSlot(void) const;
    int getNearestFaceEdgeVertexIndexA(void) const;
    int getNearestFaceEdgeVertexIndexB(void) const;

    void setNearestFaceVertex(int vertexSlot, int vertexIndex);
    int getNearestFaceVertexSlot(void) const;
    int getNearestFaceVertexIndex(void) const;

    void setModelPoint(const SbVec3f &point);
    const SbVec3f &getModelPoint(void) const;

private:
    SbString dbpath;
    SbString sourceInstanceKey;
    SbString sourceName;
    SbString sourceType;
    SbString materialShader;
    SbString editIntentId;
    SbString editIntentRole;
    SbVec3f modelPoint;
    SbColor materialColor;
    uint32_t sourceId;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    PrimitiveKind primitiveKind;
    int primitiveIndex;
    int faceVertexIndex[3];
    int nearestFaceEdgeSlot;
    int nearestFaceEdgeVertexIndex[2];
    int nearestFaceVertexSlot;
    int nearestFaceVertexIndex;
};

struct BOBOL_EXPORT BObolSourceMeshPickResult {
    SoBRLPickDetail detail;
    SbVec3f point;
    float distance;
    SbBool hit;

    BObolSourceMeshPickResult(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolRtPickResult {
    SoBRLPickDetail detail;
    SbVec3f point;
    SbVec3f normal;
    float distance;
    SbBool hit;
    /** Parametric surface coordinates; synthetic clip hits have no UV. */
    SbVec2f uv;
    SbBool uvValid;

    BObolRtPickResult(void);
    void clear(void);
};

class BOBOL_EXPORT BObolRtPickCache {
public:
    BObolRtPickCache(void);
    ~BObolRtPickCache(void);

    void clear(void);
    SbBool prepare(struct db_i *dbip,
	    const std::vector<SbString> &objectPaths);
    SbBool isReady(void) const;
    int getObjectPathCount(void) const;
    SbBool pickRay(BObolRtPickResult &pick,
	    const SbVec3f &rayOrigin,
	    const SbVec3f &rayDirection,
	    const SbPlane *clipPlanes = NULL,
	    size_t clipPlaneCount = 0,
	    struct resource *resource = NULL) const;
    /** Return the nearest ray boundary at or beyond @p minimumDistance.
     * A positive minimum lets worker-side secondary rays starting just inside
     * a solid select its exit boundary instead of repeating the entry point. */
    SbBool pickRay(BObolRtPickResult &pick,
	    const SbVec3f &rayOrigin,
	    const SbVec3f &rayDirection,
	    const SbPlane *clipPlanes,
	    size_t clipPlaneCount,
	    struct resource *resource,
	    float minimumDistance) const;

    /**
     * Initialize one librt resource for this retained acceleration structure.
     * The caller owns @p resource and must keep it alive until this cache is
     * cleared, because librt records the resource pointer in the rt_i.
     */
    SbBool initializeResource(struct resource *resource, int cpu) const;

private:
    BObolRtPickCache(const BObolRtPickCache &);
    BObolRtPickCache &operator=(const BObolRtPickCache &);

    struct rt_i *rtip;
    struct db_i *database;
    std::vector<SbString> objectPaths;
    SbBool ready;
};

BOBOL_EXPORT SbBool
bobol_pick_source_full_detail_result(
	BObolSourceMeshPickResult &pick,
	const BObolSourceMeshRequest &sourceRequest,
	const BObolLodResult &result,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	const SbPlane *clipPlanes = NULL,
	size_t clipPlaneCount = 0);

BOBOL_EXPORT SbBool
bobol_pick_rt_ray(BObolRtPickResult &pick,
	struct db_i *dbip,
	const std::vector<SbString> &objectPaths,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection);

class BOBOL_EXPORT SoBRLSourceMeshPickAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLSourceMeshPickAction);

public:
    SoBRLSourceMeshPickAction(void);
    virtual ~SoBRLSourceMeshPickAction(void);
    static void initClass(void);

    void setRay(const SbVec3f &origin, const SbVec3f &direction);
    const SbVec3f &getRayOrigin(void) const;
    const SbVec3f &getRayDirection(void) const;
    int getVisitedMeshCount(void) const;
    int getSourceBackedFullDetailRequestCount(void) const;
    const BObolSourceMeshRequest &getSourceBackedFullDetailRequest(int index) const;
    SbBool makeSourceBackedFullDetailLodRequest(int index,
	    BObolLodRequest &request,
	    const BObolLodRequest *templateRequest = 0) const;
    int submitSourceBackedFullDetailRequests(BObolLodService *service,
	    uint64_t generation, struct db_i *dbip,
	    const BObolLodRequest *templateRequest = 0,
	    uint64_t maxFullDetailFaceCount = 0,
	    uint64_t maxFullDetailPointCount = 0) const;
    int consumeSourceBackedFullDetailResults(
	    BObolSourceMeshPickResult &pick,
	    const std::vector<BObolLodResult> &results,
	    const BObolLodRequest *templateRequest = 0) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void resetResults(void);
    void appendSourceBackedFullDetailRequest(const SoBRLMeshShape *shape,
	    const SbMatrix &localToWorld);
    SbBool rayIntersectsBounds(const BObolSourceMeshRequest &request) const;

    std::vector<BObolSourceMeshRequest> sourceBackedFullDetailRequests;
    SbVec3f rayOrigin;
    SbVec3f rayDirection;
    int visitedMeshCount;
};

#endif /* BOBOL_BPICKDETAIL_H */
