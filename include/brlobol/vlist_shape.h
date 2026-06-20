/*                   V L I S T _ S H A P E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/vlist_shape.h */

#ifndef BRLOBOL_VLIST_SHAPE_H
#define BRLOBOL_VLIST_SHAPE_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoMFInt32.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoShape.h>

class SoPickedPoint;
class SoGLRenderAction;
class SoRayPickAction;

class BRLOBOL_EXPORT SoBRLVListShape : public SoShape {
    typedef SoShape inherited;

    SO_NODE_HEADER(SoBRLVListShape);

public:
    enum Command {
	MOVE = 0,
	DRAW = 1,
	POINT = 2
    };

    SoMFVec3f point;
    SoMFInt32 command;
    SoSFString sourcePath;
    SoSFString sourceName;
    SoSFString sourceType;
    SoSFUInt32 sourceId;
    SoSFInt32 regionId;
    SoSFInt32 airCode;
    SoSFInt32 materialId;
    SoSFInt32 los;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFString materialShader;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFColor selectedColor;
    SoSFColor highlightedColor;
    SoSFColor ghostedColor;
    SoSFBool visible;
    SoSFBool selectable;
    SoSFBool selected;
    SoSFBool highlighted;
    SoSFBool ghosted;
    SoSFBool hiddenLine;
    SoSFBool editEmphasis;
    SoSFUInt32 lodPolicy;
    SoMFInt32 selectedPrimitive;
    SoMFInt32 highlightedPrimitive;

    SoBRLVListShape(void);
    static void initClass(void);

    void setLineSet(const SbVec3f *points, const int32_t *commands, int count);
    int getSegmentCount(void) const;
    SbBool getSegment(int segmentIndex, SbVec3f &a, SbVec3f &b) const;
    int getPointPrimitiveCount(void) const;
    SbBool getPointPrimitive(int pointIndex, int &primitiveIndex, SbVec3f &pointOut) const;
    SbBool isPrimitiveSelected(int primitiveIndex) const;
    SbBool isPrimitiveHighlighted(int primitiveIndex) const;

    void GLRender(SoGLRenderAction *action) override;
    void computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center) override;

protected:
    virtual ~SoBRLVListShape(void);
    void generatePrimitives(SoAction *action) override;
    SoDetail *createLineSegmentDetail(SoRayPickAction *action,
	    const SoPrimitiveVertex *v1,
	    const SoPrimitiveVertex *v2,
	    SoPickedPoint *pp) override;
    SoDetail *createPointDetail(SoRayPickAction *action,
	    const SoPrimitiveVertex *v,
	    SoPickedPoint *pp) override;
};

#endif /* BRLOBOL_VLIST_SHAPE_H */
