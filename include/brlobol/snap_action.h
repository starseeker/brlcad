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

#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>

#include <stdint.h>

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
	ALL_KINDS = ENDPOINT | MIDPOINT | LINE_NEAREST | FACE_NEAREST | CENTER | CONSTRUCTION_PLANE
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
    void setConstructionPlane(const SbVec3f &origin, const SbVec3f &normal);
    void setConstructionPlane(const SbVec3f &origin, const SbVec3f &normal,
	    const SbString &path);
    void clearConstructionPlane(void);
    SbBool hasConstructionPlane(void) const;

    SbBool hasCandidate(void) const;
    SnapKind getKind(void) const;
    const SbVec3f &getPoint(void) const;
    const SbString &getPath(void) const;
    int getPrimitiveIndex(void) const;
    float getDistance(void) const;

protected:
    virtual void beginTraversal(SoNode *node);

private:
    static void nodeAction(SoAction *action, SoNode *node);
    static void vlistShapeAction(SoAction *action, SoNode *node);
    static void meshShapeAction(SoAction *action, SoNode *node);

    void consider(SnapKind kind, const SbString &path, int primitiveIndex,
	    const SbVec3f &query, const SbVec3f &candidate);
    SbVec3f pointForCoordinateSpace(const SbMatrix &localToWorld,
	    const SbVec3f &localPoint) const;
    SbBool selectionAllows(SbBool selected) const;

    SbVec3f queryPoint;
    SbVec3f candidatePoint;
    SbString candidatePath;
    uint32_t enabledKinds;
    float tolerance;
    float bestDistance;
    int candidatePrimitiveIndex;
    SnapKind candidateKind;
    SbVec3f constructionPlaneOrigin;
    SbVec3f constructionPlaneNormal;
    SbString constructionPlanePath;
    SelectionFilter selectionFilter;
    CoordinateSpace coordinateSpace;
    PriorityPolicy priorityPolicy;
    GeometryPolicy geometryPolicy;
    unsigned int skippedLodDisplayMeshCount;
    SbBool foundCandidate;
    SbBool constructionPlaneEnabled;
};

#endif /* BRLOBOL_SNAP_ACTION_H */
