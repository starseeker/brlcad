/*                    S N A P _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/snap_action.h */

#ifndef BRLOBOL_SNAP_ACTION_H
#define BRLOBOL_SNAP_ACTION_H

#include "brlobol/defines.h"
#include "brlobol/lod_realization.h"
#include "brlobol/source_mesh_request.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stdint.h>
#include <vector>

class BRLObolLodService;
class SoBRLDatabaseSource;
class SoBRLMeshShape;
struct db_i;

class BRLOBOL_EXPORT SoBRLSnapAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLSnapAction);

public:
    enum SnapKind {
	NONE = 0,
	ENDPOINT = 1,
	MIDPOINT = 2,
	LINE_NEAREST = 4,
	FACE_NEAREST = 8,
	CENTER = 16,
	CONSTRUCTION_PLANE = 32,
	VERTEX = 64,
	EDGE_NEAREST = 128,
	ALL_KINDS = ENDPOINT | MIDPOINT | LINE_NEAREST | FACE_NEAREST | CENTER | CONSTRUCTION_PLANE | VERTEX | EDGE_NEAREST
    };

    enum SelectionFilter {
	ALL_GEOMETRY = 0,
	SELECTED_ONLY = 1,
	UNSELECTED_ONLY = 2
    };

    enum CoordinateSpace {
	WORLD_SPACE = 0,
	PATH_LOCAL_SPACE = 1
    };

    enum PriorityPolicy {
	NEAREST_DISTANCE = 0,
	FEATURE_PRIORITY = 1
    };

    enum GeometryPolicy {
	FULL_DETAIL = 0,
	DISPLAY_LEVEL = 1
    };

    SoBRLSnapAction(void);
    virtual ~SoBRLSnapAction(void);
    static void initClass(void);

    void setQueryPoint(const SbVec3f &point);
    void setTolerance(float tolerance);
    void setEnabledKinds(uint32_t kinds);
    uint32_t getEnabledKinds(void) const;
    void setSelectionFilter(SelectionFilter filter);
    SelectionFilter getSelectionFilter(void) const;
    void setCoordinateSpace(CoordinateSpace space);
    CoordinateSpace getCoordinateSpace(void) const;
    void setPriorityPolicy(PriorityPolicy policy);
    PriorityPolicy getPriorityPolicy(void) const;
    void setGeometryPolicy(GeometryPolicy policy);
    GeometryPolicy getGeometryPolicy(void) const;
    unsigned int getSkippedLodDisplayMeshCount(void) const;
    int getSourceBackedFullDetailRequestCount(void) const;
    const BRLObolSourceMeshRequest &getSourceBackedFullDetailRequest(int index) const;
    SbBool makeSourceBackedFullDetailLodRequest(int index,
	    BRLObolLodRequest &request,
	    const BRLObolLodRequest *templateRequest = 0) const;
    SbBool consumeSourceBackedFullDetailResult(
	    const BRLObolSourceMeshRequest &sourceRequest,
	    const BRLObolLodResult &result);
    int submitSourceBackedFullDetailRequests(BRLObolLodService *service,
	    uint64_t generation, struct db_i *dbip,
	    const BRLObolLodRequest *templateRequest = 0,
	    uint64_t maxFullDetailFaceCount = 0,
	    uint64_t maxFullDetailPointCount = 0) const;
    int consumeSourceBackedFullDetailResults(
	    const std::vector<BRLObolLodResult> &results,
	    const BRLObolLodRequest *templateRequest = 0);
    void setConstructionPlane(const SbVec3f &origin, const SbVec3f &normal);
    void setConstructionPlane(const SbVec3f &origin, const SbVec3f &normal,
	    const SbString &path);
    void clearConstructionPlane(void);
    SbBool hasConstructionPlane(void) const;

    SbBool hasCandidate(void) const;
    SnapKind getKind(void) const;
    const SbVec3f &getPoint(void) const;
    const SbString &getPath(void) const;
    const SbString &getEditIntentId(void) const;
    const SbString &getEditIntentRole(void) const;
    int getPrimitiveIndex(void) const;
    int getVertexIndex(void) const;
    int getEdgeSlot(void) const;
    int getEdgeVertexIndexA(void) const;
    int getEdgeVertexIndexB(void) const;
    float getDistance(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    friend class SoBRLDatabaseSource;
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    static void vlistShapeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void appendSourceBackedFullDetailRequest(const SoBRLMeshShape *shape,
	    const SbMatrix &localToWorld);
    void consider(SnapKind kind, const SbString &path,
	    const SbString &editIntentId,
	    const SbString &editIntentRole,
	    int primitiveIndex,
	    const SbVec3f &query, const SbVec3f &candidate,
	    int vertexIndex = -1,
	    int edgeSlot = -1,
	    int edgeVertexIndexA = -1,
	    int edgeVertexIndexB = -1);
    SbVec3f pointForCoordinateSpace(const SbMatrix &localToWorld,
	    const SbVec3f &localPoint) const;
    SbBool selectionAllows(SbBool selected) const;

    SbVec3f queryPoint;
    SbVec3f candidatePoint;
    SbString candidatePath;
    SbString candidateEditIntentId;
    SbString candidateEditIntentRole;
    uint32_t enabledKinds;
    float tolerance;
    float bestDistance;
    int candidatePrimitiveIndex;
    int candidateVertexIndex;
    int candidateEdgeSlot;
    int candidateEdgeVertexIndex[2];
    SnapKind candidateKind;
    SbVec3f constructionPlaneOrigin;
    SbVec3f constructionPlaneNormal;
    SbString constructionPlanePath;
    SelectionFilter selectionFilter;
    CoordinateSpace coordinateSpace;
    PriorityPolicy priorityPolicy;
    GeometryPolicy geometryPolicy;
    unsigned int skippedLodDisplayMeshCount;
    std::vector<BRLObolSourceMeshRequest> sourceBackedFullDetailRequests;
    SbBool foundCandidate;
    SbBool constructionPlaneEnabled;
};

#endif /* BRLOBOL_SNAP_ACTION_H */
