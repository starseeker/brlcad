/*                  M E A S U R E _ A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/measure_action.h */

#ifndef BRLOBOL_MEASURE_ACTION_H
#define BRLOBOL_MEASURE_ACTION_H

#include "brlobol/defines.h"
#include "brlobol/lod_realization.h"
#include "brlobol/source_mesh_request.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <vector>

class BRLObolLodService;
class SoBRLMeshShape;
struct db_i;

class BRLOBOL_EXPORT SoBRLMeasureAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLMeasureAction);

public:
    enum PrimitiveKind {
	NONE = 0,
	LINE_SEGMENT = 1,
	FACE = 3
    };

    enum CoordinateSpace {
	WORLD_SPACE = 0,
	PATH_LOCAL_SPACE = 1
    };

    enum SelectionFilter {
	ALL_SELECTION = 0,
	SELECTED_ONLY = 1,
	UNSELECTED_ONLY = 2
    };

    enum HighlightFilter {
	ALL_HIGHLIGHT = 0,
	HIGHLIGHTED_ONLY = 1,
	UNHIGHLIGHTED_ONLY = 2
    };

    enum GeometryPolicy {
	FULL_DETAIL = 0,
	DISPLAY_LEVEL = 1
    };

    SoBRLMeasureAction(void);
    virtual ~SoBRLMeasureAction(void);
    static void initClass(void);

    void setQueryPoint(const SbVec3f &point);
    void clearQueryPoint(void);
    void setCoordinateSpace(CoordinateSpace space);
    CoordinateSpace getCoordinateSpace(void) const;
    void setSelectionFilter(SelectionFilter filter);
    SelectionFilter getSelectionFilter(void) const;
    void setHighlightFilter(HighlightFilter filter);
    HighlightFilter getHighlightFilter(void) const;
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

    SbBool hasSegments(void) const;
    int getShapeCount(void) const;
    int getSegmentCount(void) const;
    SbBool hasFaces(void) const;
    int getTriangleCount(void) const;
    float getSurfaceArea(void) const;
    float getTotalLength(void) const;
    const SbBox3f &getBounds(void) const;

    SbBool hasNearestSegment(void) const;
    SbBool hasNearestPrimitive(void) const;
    PrimitiveKind getNearestPrimitiveKind(void) const;
    const SbString &getNearestPath(void) const;
    int getNearestPrimitiveIndex(void) const;
    const SbVec3f &getNearestPoint(void) const;
    float getNearestDistance(void) const;
    SbBool hasAngle(void) const;
    float getAngleDegrees(void) const;
    const SbVec3f &getAnglePoint(void) const;
    const SbString &getAnglePath(void) const;
    int getAnglePrimitiveIndexA(void) const;
    int getAnglePrimitiveIndexB(void) const;
    float getAngleDistance(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void vlistShapeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void resetResults(void);
    void appendSourceBackedFullDetailRequest(const SoBRLMeshShape *shape,
	    const SbMatrix &localToWorld);
    void measureSegment(const SbString &path, int primitiveIndex,
	    const SbVec3f &a, const SbVec3f &b);
    void measureTriangle(const SbString &path, int primitiveIndex,
	    const SbVec3f &a, const SbVec3f &b, const SbVec3f &c);
    void considerAngle(const SbString &path,
	    int primitiveIndexA,
	    int primitiveIndexB,
	    const SbVec3f &point,
	    float angleDegrees);
    SbVec3f pointForCoordinateSpace(const SbMatrix &localToWorld,
	    const SbVec3f &localPoint) const;
    SbBool selectionAllows(SbBool selected) const;
    SbBool highlightAllows(SbBool highlighted) const;

    SbVec3f queryPoint;
    SbVec3f nearestPoint;
    SbVec3f anglePoint;
    SbString nearestPath;
    SbString anglePath;
    SbBox3f bounds;
    float totalLength;
    float surfaceArea;
    float nearestDistance;
    float angleDegrees;
    float angleDistance;
    int shapeCount;
    int segmentCount;
    int triangleCount;
    int nearestPrimitiveIndex;
    int anglePrimitiveIndexA;
    int anglePrimitiveIndexB;
    PrimitiveKind nearestPrimitiveKind;
    CoordinateSpace coordinateSpace;
    SelectionFilter selectionFilter;
    HighlightFilter highlightFilter;
    GeometryPolicy geometryPolicy;
    unsigned int skippedLodDisplayMeshCount;
    std::vector<BRLObolSourceMeshRequest> sourceBackedFullDetailRequests;
    SbBool haveQueryPoint;
    SbBool haveNearestPrimitive;
    SbBool haveAngle;
};

#endif /* BRLOBOL_MEASURE_ACTION_H */
