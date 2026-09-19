/*              B E D I T M A N I P U L A T O R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BEditManipulator.h */

#ifndef BOBOL_BEDITMANIPULATOR_H
#define BOBOL_BEDITMANIPULATOR_H

#include "BObol/BDefines.h"

#include <cstdint>

#include <Inventor/SbVec3f.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

class SoCamera;

/**
 * Retained, model-space handles for descriptor-driven primitive editing.
 *
 * The node owns presentation and screen-space hit/drag projection only.  It
 * deliberately has no database or transaction knowledge: a client maps its
 * stable Handle values to semantic librt edit operations and applies those
 * operations through the client's shared edit-session authority.
 */
class BOBOL_EXPORT SoBRLEditManipulator : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLEditManipulator);

public:
    enum Handle {
	HANDLE_NONE = -1,
	HANDLE_AXIS_A = 0,
	HANDLE_AXIS_B = 1,
	HANDLE_AXIS_C = 2
    };

    SoSFString manipulatorId;
    SoSFUInt32 sessionRevision;
    SoSFBool visible;
    SoSFInt32 hoverHandle;
    SoSFInt32 activeHandle;

    SoBRLEditManipulator(void);
    static void initClass(void);

    /** Replace the retained center and three model-space axis handles. */
    void setEllipsoidAxes(const SbVec3f &center, const SbVec3f &axisA,
	    const SbVec3f &axisB, const SbVec3f &axisC);

    void setVisible(SbBool value);
    void setHoverHandle(Handle handle);
    void setActiveHandle(Handle handle);

    /** Return the closest endpoint within @p radiusPixels. */
    Handle hitTest(int x, int y, int width, int height,
	    const SoCamera *camera, float radiusPixels = 12.0f) const;

    /**
     * Project a pointer onto a handle's screen-space center-to-end line.
     * @p factor is one at the current endpoint, zero at the center, and may
     * exceed one while dragging outward.
     */
    SbBool projectedScale(Handle handle, int x, int y, int width, int height,
	    const SoCamera *camera, float &factor) const;

    /** Project a point along @p handle to upper-left-origin pixel space. */
    SbBool screenPosition(Handle handle, float factor, int width, int height,
	    const SoCamera *camera, int &x, int &y) const;

    SbVec3f center(void) const;
    SbVec3f axis(Handle handle) const;

    void rebuildGeometry(void);

protected:
    virtual ~SoBRLEditManipulator(void);

private:
    SbBool project(const SbVec3f &point, int width, int height,
	    const SoCamera *camera, SbVec3f &pixel) const;

    SbVec3f editCenter;
    SbVec3f editAxes[3];
};

/**
 * Retained, indexed topology presenter for primitive editing.
 *
 * One node represents an arbitrary point/edge/face set.  Feature indices are
 * stable positions in the arrays supplied to setTopology; a client maps them
 * to semantic librt edit indices.  Like SoBRLEditManipulator, this class owns
 * only presentation and screen-space hit testing and has no database/session
 * authority.
 */
class BOBOL_EXPORT SoBRLIndexedEditManipulator : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLIndexedEditManipulator);

public:
    enum Domain {
	DOMAIN_NONE = 0,
	DOMAIN_VERTEX = 1,
	DOMAIN_EDGE = 2,
	DOMAIN_FACE = 3
    };

    SoSFString manipulatorId;
    SoSFUInt32 sessionRevision;
    SoSFBool visible;
    SoSFInt32 selectionDomain;
    SoSFInt32 selectedIndex;
    SoSFInt32 hoverIndex;
    SoSFInt32 activeIndex;

    SoBRLIndexedEditManipulator(void);
    static void initClass(void);

    /**
     * Atomically replace retained topology.  Edge indices are packed pairs;
     * face indices are packed consecutively according to faceVertexCounts.
     * When edgeFeatureIndices is non-null it contains one semantic feature
     * index per edge pair.  Multiple line pairs may therefore present and hit
     * as one curved or piecewise-linear feature.  Feature indices must be a
     * dense zero-based range.  A null array maps each edge pair one-to-one.
     * vertexFeatureCount may limit vertex picking/rendering to a prefix of the
     * point array while later points remain available to edge/face topology.
     */
    void setTopology(const SbVec3f *points, int pointCount,
	    const int32_t *edgeIndices, int edgeCount,
	    const int32_t *faceIndices, const int32_t *faceVertexCounts,
	    int faceCount, const int32_t *edgeFeatureIndices = nullptr,
	    int vertexFeatureCount = -1);

    void setVisible(SbBool value);
    void setSelectionDomain(Domain domain);
    void setSelectedIndex(int index);
    void setHoverIndex(int index);
    void setActiveIndex(int index);

    int pointCount(void) const;
    int edgeCount(void) const;
    int faceCount(void) const;

    /** Return the closest feature in @p domain within the pixel tolerance. */
    int hitTest(Domain domain, int x, int y, int width, int height,
	    const SoCamera *camera, float radiusPixels = 12.0f) const;

    /** Return a feature's model-space point/midpoint/centroid. */
    SbBool featurePosition(Domain domain, int index, SbVec3f &position) const;

    /** Project a feature's representative position to upper-left pixels. */
    SbBool screenPosition(Domain domain, int index, int width, int height,
	    const SoCamera *camera, int &x, int &y) const;

    void rebuildGeometry(void);

protected:
    virtual ~SoBRLIndexedEditManipulator(void);

private:
    class Private;
    Private *d;

    SbBool project(const SbVec3f &point, int width, int height,
	    const SoCamera *camera, SbVec3f &pixel) const;
};

#endif /* BOBOL_BEDITMANIPULATOR_H */
