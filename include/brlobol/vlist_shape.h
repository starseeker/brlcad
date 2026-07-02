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

#include <Inventor/fields/SoMFBool.h>
#include <Inventor/fields/SoMFColor.h>
#include <Inventor/fields/SoMFFloat.h>
#include <Inventor/fields/SoMFString.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoMFInt32.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoShape.h>

#include <vector>

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

    enum AnnotationSegmentKind {
	ANNOTATION_SEGMENT_NONE = 0,
	ANNOTATION_SEGMENT_LINE = 1,
	ANNOTATION_SEGMENT_TEXT = 2
    };

    SoMFVec3f point;
    SoMFInt32 command;
    SoSFVec3f annotationBasePoint;
    SoMFVec3f annotationPoint;
    SoMFBool annotationSegmentTextValid;
    SoMFInt32 annotationSegmentKind;
    SoMFInt32 annotationSegmentStart;
    SoMFInt32 annotationSegmentEnd;
    SoMFInt32 annotationTextRefPoint;
    SoMFString annotationText;
    SoMFBool pointColorValid;
    SoMFColor pointColor;
    SoMFBool pointScaleValid;
    SoMFFloat pointScale;
    SoMFBool pointNormalValid;
    SoMFVec3f pointNormal;
    SoSFString sourcePath;
    SoSFString sourceName;
    SoSFString sourceType;
    SoSFUInt32 sourceId;
    SoSFString displayName;
    SoSFString geometryName;
    SoSFString cacheIdentity;
    SoSFString sourceIdentity;
    SoSFString ownerSourcePath;
    SoSFString ownerSourceInstanceKey;
    SoSFUInt32 ownerSourceRevision;
    SoSFUInt32 ownerInputsRevision;
    SoSFUInt32 ownerViewRevision;
    SoSFUInt32 ownerRealizedRevision;
    SoSFUInt32 ownerRealizedSourceRevision;
    SoSFUInt32 ownerRealizedInputsRevision;
    SoSFUInt32 ownerRealizedViewRevision;
    SoSFInt32 ownerRealizationStatus;
    SoSFString ownerRealizationDiagnostic;
    SoSFString ownerRealizationIdentity;
    SoSFBool ownerSourceStale;
    SoSFUInt32 ownerStaleReason;
    SoSFBool databaseIntent;
    SoSFBool overlayIntent;
    SoSFBool hudIntent;
    SoSFBool localSource;
    SoSFBool sharedSource;
    SoSFBool nonDatabaseSource;
    SoSFInt32 drawMode;
    SoSFString recordRole;
    SoSFString geometryKind;
    SoSFInt32 regionId;
    SoSFInt32 airCode;
    SoSFInt32 materialId;
    SoSFInt32 los;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFString materialShader;
    SoSFUInt32 materialRevision;
    SoSFBool drawMatrixValid;
    SoSFMatrix drawMatrix;
    SoSFBool drawCenterValid;
    SoSFVec3f drawCenter;
    SoSFBool drawSizeValid;
    SoSFFloat drawSize;
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
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool hiddenLine;
    SoSFBool editEmphasis;
    SoSFString editIntentId;
    SoSFString editIntentRole;
    SoSFUInt32 lodPolicy;
    SoMFInt32 selectedPrimitive;
    SoMFInt32 highlightedPrimitive;

    SoBRLVListShape(void);
    static void initClass(void);

    void setLineSet(const SbVec3f *points, const int32_t *commands, int count);
    void setPrecisePoints(const double *points, int count);
    SbBool getPrecisePoint(int index, double *pointOut) const;
    void setPreciseAnnotationPoints(const double *points, int count);
    SbBool getPreciseAnnotationPoint(int index, double *pointOut) const;
    SbBool translatePoints(const SbVec3f &offset);
    void setDrawCenter(const SbVec3f &center);
    SbBool updateDrawBoundsFromPoints(void);
    void setPointAttributes(const int *colorValid, const SbColor *colors,
	    const int *scaleValid, const float *scales,
	    const int *normalValid, const SbVec3f *normals, int count);
    int getSegmentCount(void) const;
    SbBool getSegment(int segmentIndex, SbVec3f &a, SbVec3f &b) const;
    int getPointPrimitiveCount(void) const;
    SbBool getPointPrimitive(int pointIndex, int &primitiveIndex, SbVec3f &pointOut) const;
    SbBool getPointColor(int primitiveIndex, SbColor &colorOut) const;
    SbBool getPointScale(int primitiveIndex, float &scaleOut) const;
    SbBool getPointNormal(int primitiveIndex, SbVec3f &normalOut) const;
    SbBool isPrimitiveSelected(int primitiveIndex) const;
    SbBool isPrimitiveHighlighted(int primitiveIndex) const;

    void GLRender(SoGLRenderAction *action) override;
    void computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center) override;

private:
    std::vector<double> precisePoints;
    std::vector<double> preciseAnnotationPoints;

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
