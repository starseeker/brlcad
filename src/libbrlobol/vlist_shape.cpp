/*                   V L I S T _ S H A P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_realization.h"
#include "brlobol/pick_detail.h"
#include "brlobol/vlist_shape.h"

#include <obol/cad/SoCADViewState.h>

#include <Inventor/SoDB.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoPointDetail.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLDisplayList.h>
#include <Inventor/elements/SoContextManagerElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoProjectionMatrixElement.h>
#include <Inventor/elements/SoViewingMatrixElement.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/gl.h>

#include <algorithm>
#include <cmath>
#include <vector>

SO_NODE_SOURCE(SoBRLVListShape);

SoBRLVListShape::SoBRLVListShape(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLVListShape);

    SO_NODE_DEFINE_ENUM_VALUE(Command, MOVE);
    SO_NODE_DEFINE_ENUM_VALUE(Command, DRAW);
    SO_NODE_DEFINE_ENUM_VALUE(Command, POINT);

    SO_NODE_ADD_EMPTY_MFIELD(point);
    SO_NODE_ADD_EMPTY_MFIELD(command);
    SO_NODE_ADD_FIELD(annotationBasePoint, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_EMPTY_MFIELD(annotationPoint);
    SO_NODE_ADD_EMPTY_MFIELD(annotationSegmentTextValid);
    SO_NODE_ADD_EMPTY_MFIELD(annotationSegmentKind);
    SO_NODE_ADD_EMPTY_MFIELD(annotationSegmentStart);
    SO_NODE_ADD_EMPTY_MFIELD(annotationSegmentEnd);
    SO_NODE_ADD_EMPTY_MFIELD(annotationTextRefPoint);
    SO_NODE_ADD_EMPTY_MFIELD(annotationText);
    SO_NODE_ADD_EMPTY_MFIELD(pointColorValid);
    SO_NODE_ADD_EMPTY_MFIELD(pointColor);
    SO_NODE_ADD_EMPTY_MFIELD(pointScaleValid);
    SO_NODE_ADD_EMPTY_MFIELD(pointScale);
    SO_NODE_ADD_EMPTY_MFIELD(pointNormalValid);
    SO_NODE_ADD_EMPTY_MFIELD(pointNormal);
    SO_NODE_ADD_FIELD(sourcePath, (""));
    SO_NODE_ADD_FIELD(sourceName, (""));
    SO_NODE_ADD_FIELD(sourceType, (""));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(displayName, (""));
    SO_NODE_ADD_FIELD(geometryName, (""));
    SO_NODE_ADD_FIELD(cacheIdentity, (""));
    SO_NODE_ADD_FIELD(sourceIdentity, (""));
    SO_NODE_ADD_FIELD(ownerSourcePath, (""));
    SO_NODE_ADD_FIELD(ownerSourceInstanceKey, (""));
    SO_NODE_ADD_FIELD(ownerSourceRevision, (0));
    SO_NODE_ADD_FIELD(ownerInputsRevision, (0));
    SO_NODE_ADD_FIELD(ownerViewRevision, (0));
    SO_NODE_ADD_FIELD(ownerRealizedRevision, (0));
    SO_NODE_ADD_FIELD(ownerRealizedSourceRevision, (0));
    SO_NODE_ADD_FIELD(ownerRealizedInputsRevision, (0));
    SO_NODE_ADD_FIELD(ownerRealizedViewRevision, (0));
    SO_NODE_ADD_FIELD(ownerRealizationStatus, (0));
    SO_NODE_ADD_FIELD(ownerRealizationDiagnostic, (""));
    SO_NODE_ADD_FIELD(ownerRealizationIdentity, (""));
    SO_NODE_ADD_FIELD(ownerSourceStale, (FALSE));
    SO_NODE_ADD_FIELD(ownerStaleReason, (0));
    SO_NODE_ADD_FIELD(databaseIntent, (FALSE));
    SO_NODE_ADD_FIELD(overlayIntent, (FALSE));
    SO_NODE_ADD_FIELD(hudIntent, (FALSE));
    SO_NODE_ADD_FIELD(localSource, (FALSE));
    SO_NODE_ADD_FIELD(sharedSource, (FALSE));
    SO_NODE_ADD_FIELD(nonDatabaseSource, (FALSE));
    SO_NODE_ADD_FIELD(drawMode, (0));
    SO_NODE_ADD_FIELD(recordRole, (""));
    SO_NODE_ADD_FIELD(geometryKind, (""));
    SO_NODE_ADD_FIELD(regionId, (0));
    SO_NODE_ADD_FIELD(airCode, (0));
    SO_NODE_ADD_FIELD(materialId, (0));
    SO_NODE_ADD_FIELD(los, (0));
    SO_NODE_ADD_FIELD(materialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(materialColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialShader, (""));
    SO_NODE_ADD_FIELD(materialRevision, (0));
    SO_NODE_ADD_FIELD(drawMatrixValid, (FALSE));
    SO_NODE_ADD_FIELD(drawMatrix, (SbMatrix::identity()));
    SO_NODE_ADD_FIELD(drawCenterValid, (FALSE));
    SO_NODE_ADD_FIELD(drawCenter, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(drawSizeValid, (FALSE));
    SO_NODE_ADD_FIELD(drawSize, (0.0f));
    SO_NODE_ADD_FIELD(colorOverride, (FALSE));
    SO_NODE_ADD_FIELD(color, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(selectedColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(highlightedColor, (SbColor(1.0f, 1.0f, 0.0f)));
    SO_NODE_ADD_FIELD(ghostedColor, (SbColor(0.55f, 0.55f, 0.55f)));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selectable, (TRUE));
    SO_NODE_ADD_FIELD(selected, (FALSE));
    SO_NODE_ADD_FIELD(highlighted, (FALSE));
    SO_NODE_ADD_FIELD(ghosted, (FALSE));
    SO_NODE_ADD_FIELD(lineStyle, (0));
    SO_NODE_ADD_FIELD(lineWidth, (0));
    SO_NODE_ADD_FIELD(transparency, (0.0f));
    SO_NODE_ADD_FIELD(hiddenLine, (FALSE));
    SO_NODE_ADD_FIELD(editEmphasis, (FALSE));
    SO_NODE_ADD_FIELD(editIntentId, (""));
    SO_NODE_ADD_FIELD(editIntentRole, (""));
    SO_NODE_ADD_FIELD(lodPolicy, (0));
    SO_NODE_ADD_EMPTY_MFIELD(selectedPrimitive);
    SO_NODE_ADD_EMPTY_MFIELD(highlightedPrimitive);
    SO_NODE_ADD_FIELD(sharedGeometry, (NULL));
}

SoBRLVListShape::~SoBRLVListShape(void)
{
    this->clearRenderLists();
}

void
SoBRLVListShape::clearRenderLists(void)
{
    for (std::map<int, SoGLDisplayList *>::iterator it =
	this->renderLists.begin(); it != this->renderLists.end(); ++it) {
	if (it->second)
	    it->second->unref();
    }
    this->renderLists.clear();
}

void
SoBRLVListShape::notify(SoNotList *list)
{
    this->clearRenderLists();
    inherited::notify(list);
}

void
SoBRLVListShape::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLVListShape, SoShape, "Shape");
}

void
SoBRLVListShape::setSharedGeometry(SoBRLVListShape *shape)
{
    this->sharedGeometry = shape;
}

SoBRLVListShape *
SoBRLVListShape::getSharedGeometrySource(void)
{
    SoNode *node = this->sharedGeometry.getValue();
    if (node && node != this &&
	node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<SoBRLVListShape *>(node);
    return this;
}

const SoBRLVListShape *
SoBRLVListShape::getSharedGeometrySource(void) const
{
    const SoNode *node = this->sharedGeometry.getValue();
    if (node && node != this &&
	node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node);
    return this;
}

SoBRLVListShape *
SoBRLVListShape::getGeometrySource(void)
{
    return this->getSharedGeometrySource();
}

const SoBRLVListShape *
SoBRLVListShape::getGeometrySource(void) const
{
    return this->getSharedGeometrySource();
}

void
SoBRLVListShape::setLineSet(const SbVec3f *points, const int32_t *commands, int count)
{
    this->sharedGeometry = NULL;
    this->precisePoints.clear();
    this->preciseAnnotationPoints.clear();
    this->point.setNum(0);
    this->command.setNum(0);
    this->annotationPoint.setNum(0);
    this->annotationSegmentTextValid.setNum(0);
    this->annotationSegmentKind.setNum(0);
    this->annotationSegmentStart.setNum(0);
    this->annotationSegmentEnd.setNum(0);
    this->annotationTextRefPoint.setNum(0);
    this->annotationText.setNum(0);
    this->pointColorValid.setNum(0);
    this->pointColor.setNum(0);
    this->pointScaleValid.setNum(0);
    this->pointScale.setNum(0);
    this->pointNormalValid.setNum(0);
    this->pointNormal.setNum(0);
    if (!points || !commands || count <= 0)
	return;

    this->point.setValues(0, count, points);
    this->command.setValues(0, count, commands);
}

void
SoBRLVListShape::setPrecisePoints(const double *points, int count)
{
    this->clearRenderLists();
    this->precisePoints.clear();
    if (!points || count <= 0)
	return;

    this->precisePoints.resize(static_cast<size_t>(count) * 3);
    for (int i = 0; i < count; i++) {
	this->precisePoints[static_cast<size_t>(i) * 3 + 0] =
	    points[static_cast<size_t>(i) * 3 + 0];
	this->precisePoints[static_cast<size_t>(i) * 3 + 1] =
	    points[static_cast<size_t>(i) * 3 + 1];
	this->precisePoints[static_cast<size_t>(i) * 3 + 2] =
	    points[static_cast<size_t>(i) * 3 + 2];
    }
}

SbBool
SoBRLVListShape::getPrecisePoint(int index, double *pointOut) const
{
    if (!pointOut || index < 0)
	return FALSE;

    const SoBRLVListShape *geom = this->getGeometrySource();
    const size_t offset = static_cast<size_t>(index) * 3;
    if (offset + 2 >= geom->precisePoints.size())
	return FALSE;

    pointOut[0] = geom->precisePoints[offset + 0];
    pointOut[1] = geom->precisePoints[offset + 1];
    pointOut[2] = geom->precisePoints[offset + 2];
    return TRUE;
}

void
SoBRLVListShape::setPreciseAnnotationPoints(const double *points, int count)
{
    this->preciseAnnotationPoints.clear();
    this->annotationPoint.setNum(0);
    if (!points || count <= 0)
	return;

    this->preciseAnnotationPoints.resize(static_cast<size_t>(count) * 3);
    std::vector<SbVec3f> floatPoints(static_cast<size_t>(count));
    for (int i = 0; i < count; i++) {
	const size_t offset = static_cast<size_t>(i) * 3;
	this->preciseAnnotationPoints[offset + 0] = points[offset + 0];
	this->preciseAnnotationPoints[offset + 1] = points[offset + 1];
	this->preciseAnnotationPoints[offset + 2] = points[offset + 2];
	floatPoints[static_cast<size_t>(i)] = SbVec3f(
		static_cast<float>(points[offset + 0]),
		static_cast<float>(points[offset + 1]),
		static_cast<float>(points[offset + 2]));
    }
    this->annotationPoint.setValues(0, count, floatPoints.data());
}

SbBool
SoBRLVListShape::getPreciseAnnotationPoint(int index, double *pointOut) const
{
    if (!pointOut || index < 0)
	return FALSE;

    const SoBRLVListShape *geom = this->getGeometrySource();
    const size_t offset = static_cast<size_t>(index) * 3;
    if (offset + 2 < geom->preciseAnnotationPoints.size()) {
	pointOut[0] = geom->preciseAnnotationPoints[offset + 0];
	pointOut[1] = geom->preciseAnnotationPoints[offset + 1];
	pointOut[2] = geom->preciseAnnotationPoints[offset + 2];
	return TRUE;
    }

    if (index >= geom->annotationPoint.getNum())
	return FALSE;

    const SbVec3f &pointValue = geom->annotationPoint[index];
    pointOut[0] = pointValue[0];
    pointOut[1] = pointValue[1];
    pointOut[2] = pointValue[2];
    return TRUE;
}

SbBool
SoBRLVListShape::translatePoints(const SbVec3f &offset)
{
    const int count = this->point.getNum();
    if (count <= 0)
	return FALSE;

    for (int i = 0; i < count; i++) {
	const SbVec3f pointValue = this->point[i];
	this->point.set1Value(i, pointValue + offset);
	const size_t preciseOffset = static_cast<size_t>(i) * 3;
	if (preciseOffset + 2 < this->precisePoints.size()) {
	    this->precisePoints[preciseOffset + 0] += offset[0];
	    this->precisePoints[preciseOffset + 1] += offset[1];
	    this->precisePoints[preciseOffset + 2] += offset[2];
	}
    }

    if (this->drawCenterValid.getValue()) {
	this->drawCenter = this->drawCenter.getValue() + offset;
    }

    return TRUE;
}

void
SoBRLVListShape::setDrawCenter(const SbVec3f &center)
{
    this->drawCenterValid = TRUE;
    this->drawCenter = center;
}

SbBool
SoBRLVListShape::updateDrawBoundsFromPoints(void)
{
    const int count = this->point.getNum();
    if (count <= 0) {
	this->drawCenterValid = FALSE;
	this->drawSizeValid = FALSE;
	return FALSE;
    }

    SbVec3f bmin = this->point[0];
    SbVec3f bmax = this->point[0];
    for (int i = 1; i < count; i++) {
	const SbVec3f pointValue = this->point[i];
	for (int axis = 0; axis < 3; axis++) {
	    if (pointValue[axis] < bmin[axis])
		bmin[axis] = pointValue[axis];
	    if (pointValue[axis] > bmax[axis])
		bmax[axis] = pointValue[axis];
	}
    }

    this->drawCenterValid = TRUE;
    this->drawCenter = SbVec3f(
			   (bmin[0] + bmax[0]) * 0.5f,
			   (bmin[1] + bmax[1]) * 0.5f,
			   (bmin[2] + bmax[2]) * 0.5f);

    float size = bmax[0] - bmin[0];
    if ((bmax[1] - bmin[1]) > size)
	size = bmax[1] - bmin[1];
    if ((bmax[2] - bmin[2]) > size)
	size = bmax[2] - bmin[2];
    this->drawSizeValid = TRUE;
    this->drawSize = size;
    return TRUE;
}

void
SoBRLVListShape::setPointAttributes(const int *colorValid,
				    const SbColor *colors, const int *scaleValid, const float *scales,
				    const int *normalValid, const SbVec3f *normals, int count)
{
    this->pointColorValid.setNum(0);
    this->pointColor.setNum(0);
    this->pointScaleValid.setNum(0);
    this->pointScale.setNum(0);
    this->pointNormalValid.setNum(0);
    this->pointNormal.setNum(0);
    if (count <= 0)
	return;

    if (colorValid && colors) {
	this->pointColorValid.setNum(count);
	for (int i = 0; i < count; i++)
	    this->pointColorValid.set1Value(i, colorValid[i] ? TRUE : FALSE);
	this->pointColor.setValues(0, count, colors);
    }
    if (scaleValid && scales) {
	this->pointScaleValid.setNum(count);
	for (int i = 0; i < count; i++)
	    this->pointScaleValid.set1Value(i, scaleValid[i] ? TRUE : FALSE);
	this->pointScale.setValues(0, count, scales);
    }
    if (normalValid && normals) {
	this->pointNormalValid.setNum(count);
	for (int i = 0; i < count; i++)
	    this->pointNormalValid.set1Value(i, normalValid[i] ? TRUE : FALSE);
	this->pointNormal.setValues(0, count, normals);
    }
}

int
SoBRLVListShape::getSegmentCount(void) const
{
    int ret = 0;
    SbVec3f last;
    SbBool haveLast = FALSE;
    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    for (int i = 0; i < n; i++) {
	switch (geom->command[i]) {
	    case MOVE:
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast)
		    ret++;
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    default:
		break;
	}
    }

    return ret;
}

SbBool
SoBRLVListShape::getSegment(int segmentIndex, SbVec3f &a, SbVec3f &b) const
{
    if (segmentIndex < 0)
	return FALSE;

    int currentSegment = 0;
    SbVec3f last;
    SbBool haveLast = FALSE;
    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    for (int i = 0; i < n; i++) {
	switch (geom->command[i]) {
	    case MOVE:
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast) {
		    if (currentSegment == segmentIndex) {
			a = last;
			b = geom->point[i];
			return TRUE;
		    }
		    currentSegment++;
		}
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    default:
		break;
	}
    }

    return FALSE;
}

int
SoBRLVListShape::getPointPrimitiveCount(void) const
{
    int ret = 0;
    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    for (int i = 0; i < n; i++) {
	if (geom->command[i] == POINT)
	    ret++;
    }

    return ret;
}

SbBool
SoBRLVListShape::getPointPrimitive(int pointIndex, int &primitiveIndex, SbVec3f &pointOut) const
{
    primitiveIndex = -1;
    pointOut = SbVec3f(0.0f, 0.0f, 0.0f);
    if (pointIndex < 0)
	return FALSE;

    int currentPoint = 0;
    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    for (int i = 0; i < n; i++) {
	if (geom->command[i] != POINT)
	    continue;
	if (currentPoint == pointIndex) {
	    primitiveIndex = i;
	    pointOut = geom->point[i];
	    return TRUE;
	}
	currentPoint++;
    }

    return FALSE;
}

SbBool
SoBRLVListShape::getPointColor(int primitiveIndex, SbColor &colorOut) const
{
    colorOut = SbColor(1.0f, 1.0f, 1.0f);
    const SoBRLVListShape *geom = this->getGeometrySource();
    if (primitiveIndex < 0 ||
	primitiveIndex >= geom->pointColorValid.getNum() ||
	primitiveIndex >= geom->pointColor.getNum() ||
	!geom->pointColorValid[primitiveIndex])
	return FALSE;

    colorOut = geom->pointColor[primitiveIndex];
    return TRUE;
}

SbBool
SoBRLVListShape::getPointScale(int primitiveIndex, float &scaleOut) const
{
    scaleOut = 0.0f;
    const SoBRLVListShape *geom = this->getGeometrySource();
    if (primitiveIndex < 0 ||
	primitiveIndex >= geom->pointScaleValid.getNum() ||
	primitiveIndex >= geom->pointScale.getNum() ||
	!geom->pointScaleValid[primitiveIndex])
	return FALSE;

    scaleOut = geom->pointScale[primitiveIndex];
    return TRUE;
}

SbBool
SoBRLVListShape::getPointNormal(int primitiveIndex, SbVec3f &normalOut) const
{
    normalOut = SbVec3f(0.0f, 0.0f, 1.0f);
    const SoBRLVListShape *geom = this->getGeometrySource();
    if (primitiveIndex < 0 ||
	primitiveIndex >= geom->pointNormalValid.getNum() ||
	primitiveIndex >= geom->pointNormal.getNum() ||
	!geom->pointNormalValid[primitiveIndex])
	return FALSE;

    normalOut = geom->pointNormal[primitiveIndex];
    return TRUE;
}

static SbBool
vlist_int_field_contains(const SoMFInt32 &field, int value)
{
    for (int i = 0; i < field.getNum(); i++) {
	if (field[i] == value)
	    return TRUE;
    }
    return FALSE;
}

SbBool
SoBRLVListShape::isPrimitiveSelected(int primitiveIndex) const
{
    return this->selected.getValue() || vlist_int_field_contains(this->selectedPrimitive, primitiveIndex);
}

SbBool
SoBRLVListShape::isPrimitiveHighlighted(int primitiveIndex) const
{
    return this->highlighted.getValue() || vlist_int_field_contains(this->highlightedPrimitive, primitiveIndex);
}

static SbBool
set_vlist_gl_state_color(SoBRLVListShape *shape, int primitiveIndex)
{
    if (shape->isPrimitiveHighlighted(primitiveIndex)) {
	const SbColor &c = shape->highlightedColor.getValue();
	glColor3f(c[0], c[1], c[2]);
	return TRUE;
    } else if (shape->isPrimitiveSelected(primitiveIndex)) {
	const SbColor &c = shape->selectedColor.getValue();
	glColor3f(c[0], c[1], c[2]);
	return TRUE;
    } else if (shape->ghosted.getValue()) {
	const SbColor &c = shape->ghostedColor.getValue();
	glColor4f(c[0], c[1], c[2], 0.35f);
	return TRUE;
    } else if (shape->colorOverride.getValue()) {
	const SbColor &c = shape->color.getValue();
	glColor3f(c[0], c[1], c[2]);
	return TRUE;
    }
    return FALSE;
}

static void
set_vlist_default_gl_color(SoBRLVListShape *shape)
{
    const SbColor &c = shape->materialColorValid.getValue() ?
		       shape->materialColor.getValue() : shape->color.getValue();
    glColor3f(c[0], c[1], c[2]);
}

static void
set_vlist_gl_color(SoBRLVListShape *shape, int primitiveIndex)
{
    if (!set_vlist_gl_state_color(shape, primitiveIndex))
	set_vlist_default_gl_color(shape);
}

static SbBool
vlist_needs_independent_segment_rendering(const SoBRLVListShape *shape)
{
    return shape &&
	   (shape->selectedPrimitive.getNum() > 0 ||
	    shape->highlightedPrimitive.getNum() > 0);
}

struct VListClipPoint {
    double v[4];
};

static VListClipPoint
vlist_clip_point(const SbMatrix &matrix, const SbVec3f &point)
{
    VListClipPoint result;
    for (int col = 0; col < 4; col++) {
	result.v[col] = static_cast<double>(point[0]) * matrix[0][col] +
	    static_cast<double>(point[1]) * matrix[1][col] +
	    static_cast<double>(point[2]) * matrix[2][col] + matrix[3][col];
    }
    return result;
}

static double
vlist_clip_plane_value(const VListClipPoint &point, int plane)
{
    switch (plane) {
	case 0: return point.v[3] + point.v[0];
	case 1: return point.v[3] - point.v[0];
	case 2: return point.v[3] + point.v[1];
	case 3: return point.v[3] - point.v[1];
	case 4: return point.v[3] + point.v[2];
	default: return point.v[3] - point.v[2];
    }
}

static SbBool
vlist_clip_segment(VListClipPoint &a, VListClipPoint &b)
{
    for (int plane = 0; plane < 6; plane++) {
	double da = vlist_clip_plane_value(a, plane);
	double db = vlist_clip_plane_value(b, plane);
	if (da < 0.0 && db < 0.0)
	    return FALSE;
	if (da >= 0.0 && db >= 0.0)
	    continue;

	const double denominator = da - db;
	if (fabs(denominator) < 1.0e-20)
	    return FALSE;
	const double t = da / denominator;
	VListClipPoint clipped;
	for (int i = 0; i < 4; i++)
	    clipped.v[i] = a.v[i] + t * (b.v[i] - a.v[i]);
	if (da < 0.0)
	    a = clipped;
	else
	    b = clipped;
    }
    return fabs(a.v[3]) > 1.0e-20 && fabs(b.v[3]) > 1.0e-20;
}

static void
vlist_put_software_pixel(unsigned char *pixels, unsigned int width,
	unsigned int height, int x, int y, const unsigned char color[4])
{
    if (!pixels || x < 0 || y < 0 ||
	static_cast<unsigned int>(x) >= width ||
	static_cast<unsigned int>(y) >= height)
	return;

    unsigned char *pixel = pixels +
	(static_cast<size_t>(y) * width + static_cast<unsigned int>(x)) * 4;
    if (color[3] == 255) {
	pixel[0] = color[0];
	pixel[1] = color[1];
	pixel[2] = color[2];
	pixel[3] = 255;
	return;
    }

    const unsigned int alpha = color[3];
    const unsigned int inverse = 255 - alpha;
    for (int i = 0; i < 3; i++)
	pixel[i] = static_cast<unsigned char>(
	    (color[i] * alpha + pixel[i] * inverse + 127) / 255);
    pixel[3] = 255;
}

static void
vlist_draw_software_line(unsigned char *pixels, unsigned int width,
	unsigned int height, int x0, int y0, int x1, int y1,
	const unsigned char color[4])
{
    int dx = abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;
    for (;;) {
	vlist_put_software_pixel(pixels, width, height, x0, y0, color);
	if (x0 == x1 && y0 == y1)
	    break;
	const int twiceError = 2 * error;
	if (twiceError >= dy) {
	    error += dy;
	    x0 += sx;
	}
	if (twiceError <= dx) {
	    error += dx;
	    y0 += sy;
	}
    }
}

static void
vlist_software_color(const SoBRLVListShape *shape, unsigned char color[4])
{
    SbColor c;
    float alpha = 1.0f;
    if (shape->highlighted.getValue())
	c = shape->highlightedColor.getValue();
    else if (shape->selected.getValue())
	c = shape->selectedColor.getValue();
    else if (shape->ghosted.getValue()) {
	c = shape->ghostedColor.getValue();
	alpha = 0.35f;
    } else if (shape->colorOverride.getValue())
	c = shape->color.getValue();
    else
	c = shape->materialColorValid.getValue() ?
	    shape->materialColor.getValue() : shape->color.getValue();

    color[0] = static_cast<unsigned char>(std::clamp(c[0], 0.0f, 1.0f) * 255.0f);
    color[1] = static_cast<unsigned char>(std::clamp(c[1], 0.0f, 1.0f) * 255.0f);
    color[2] = static_cast<unsigned char>(std::clamp(c[2], 0.0f, 1.0f) * 255.0f);
    color[3] = static_cast<unsigned char>(alpha * 255.0f);
}

static SbBool
vlist_render_software_wireframe(SoBRLVListShape *shape,
	const SoBRLVListShape *geom, int n, SoState *state)
{
    const obol::CadViewState viewState =
	SoCADViewStateElement::get(state);
    if (viewState.softwareWireMode != obol::CadSoftwareWireMode::FAST)
	return FALSE;
    if (!shape->databaseIntent.getValue() ||
	shape->drawMode.getValue() != BRLOBOL_LOD_DRAW_WIRE ||
	shape->hiddenLine.getValue() || shape->lineWidth.getValue() > 1 ||
	vlist_needs_independent_segment_rendering(shape))
	return FALSE;

    SoDB::ContextManager *manager = SoContextManagerElement::get(state);
    unsigned char *pixels = NULL;
    unsigned int width = 0, height = 0, components = 0;
    const SbBool framebuffer = manager ? manager->getCurrentSoftwareFramebuffer(
	pixels, width, height, components) : FALSE;
    if (!manager || !framebuffer || components != 4)
	return FALSE;

    const SbViewportRegion &viewport = SoViewportRegionElement::get(state);
    const SbVec2s origin = viewport.getViewportOriginPixels();
    const SbVec2s size = viewport.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return FALSE;

    const SbMatrix transform = SoModelMatrixElement::get(state) *
	SoViewingMatrixElement::get(state) *
	SoProjectionMatrixElement::get(state);
    unsigned char color[4];
    vlist_software_color(shape, color);

    VListClipPoint previous = {};
    SbBool havePrevious = FALSE;
    for (int i = 0; i < n; i++) {
	if (geom->command[i] == SoBRLVListShape::MOVE) {
	    previous = vlist_clip_point(transform, geom->point[i]);
	    havePrevious = TRUE;
	    continue;
	}
	if (geom->command[i] != SoBRLVListShape::DRAW || !havePrevious)
	    continue;

	VListClipPoint current = vlist_clip_point(transform, geom->point[i]);
	VListClipPoint a = previous;
	VListClipPoint b = current;
	previous = current;
	if (!vlist_clip_segment(a, b))
	    continue;
	const double ax = a.v[0] / a.v[3];
	const double ay = a.v[1] / a.v[3];
	const double bx = b.v[0] / b.v[3];
	const double by = b.v[1] / b.v[3];
	const int x0 = origin[0] + static_cast<int>(
	    std::lround((ax * 0.5 + 0.5) * (size[0] - 1)));
	const int y0 = origin[1] + static_cast<int>(
	    std::lround((ay * 0.5 + 0.5) * (size[1] - 1)));
	const int x1 = origin[0] + static_cast<int>(
	    std::lround((bx * 0.5 + 0.5) * (size[0] - 1)));
	const int y1 = origin[1] + static_cast<int>(
	    std::lround((by * 0.5 + 0.5) * (size[1] - 1)));
	vlist_draw_software_line(pixels, width, height, x0, y0, x1, y1, color);
    }
    return TRUE;
}

static void
vlist_gl_vertex_at(const SoBRLVListShape *shape,
		   const SoBRLVListShape *geom,
		   int index)
{
    double precisePoint[3] = {0.0, 0.0, 0.0};
    if (shape && shape->getPrecisePoint(index, precisePoint)) {
	glVertex3dv(precisePoint);
	return;
    }

    const SbVec3f &point = geom->point[index];
    glVertex3f(point[0], point[1], point[2]);
}

static void
vlist_render_independent_segments(SoBRLVListShape *shape,
				  const SoBRLVListShape *geom,
				  int n)
{
    SbBool haveLast = FALSE;
    int lastIndex = -1;
    int segmentIndex = 0;

    glBegin(GL_LINES);
    for (int i = 0; i < n; i++) {
	switch (geom->command[i]) {
	    case SoBRLVListShape::MOVE:
		lastIndex = i;
		haveLast = TRUE;
		break;
	    case SoBRLVListShape::DRAW:
		if (haveLast) {
		    set_vlist_gl_color(shape, segmentIndex);
		    vlist_gl_vertex_at(shape, geom, lastIndex);
		    vlist_gl_vertex_at(shape, geom, i);
		    segmentIndex++;
		}
		lastIndex = i;
		haveLast = TRUE;
		break;
	    default:
		break;
	}
    }
    glEnd();
}

static void
vlist_render_line_strips(SoBRLVListShape *shape,
			 const SoBRLVListShape *geom,
			 int n)
{
    SbBool stripOpen = FALSE;

    for (int i = 0; i < n; i++) {
	switch (geom->command[i]) {
	    case SoBRLVListShape::MOVE:
		if (stripOpen)
		    glEnd();
		glBegin(GL_LINE_STRIP);
		vlist_gl_vertex_at(shape, geom, i);
		stripOpen = TRUE;
		break;
	    case SoBRLVListShape::DRAW:
		if (!stripOpen) {
		    glBegin(GL_LINE_STRIP);
		    stripOpen = TRUE;
		}
		vlist_gl_vertex_at(shape, geom, i);
		break;
	    default:
		if (stripOpen) {
		    glEnd();
		    stripOpen = FALSE;
		}
		break;
	}
    }

    if (stripOpen)
	glEnd();
}

void
SoBRLVListShape::GLRender(SoGLRenderAction *action)
{
    if (!this->visible.getValue() || !this->shouldGLRender(action))
	return;

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LINE_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);
    set_vlist_gl_color(this, -1);
    if (this->lineWidth.getValue() > 0)
	glLineWidth(static_cast<GLfloat>(this->lineWidth.getValue()));

    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    if (vlist_render_software_wireframe(this, geom, n, action->getState())) {
	/* The software backend wrote ordinary database wireframes directly. */
    } else if (vlist_needs_independent_segment_rendering(this))
	vlist_render_independent_segments(this, geom, n);
    else {
	SoBRLVListShape *cacheOwner = const_cast<SoBRLVListShape *>(geom);
	SoState *state = action->getState();
	const int context = SoGLCacheContextElement::get(state);
	std::map<int, SoGLDisplayList *>::iterator found =
	    cacheOwner->renderLists.find(context);
	if (found != cacheOwner->renderLists.end()) {
	    found->second->call(state);
	} else {
	    SoGLDisplayList *list = new SoGLDisplayList(
		state, SoGLDisplayList::DISPLAY_LIST);
	    list->ref();
	    cacheOwner->renderLists[context] = list;
	    list->open(state);
	    vlist_render_line_strips(this, geom, n);
	    list->close(state);
	}
    }

    glBegin(GL_POINTS);
    for (int i = 0; i < n; i++) {
	if (geom->command[i] == POINT) {
	    if (!set_vlist_gl_state_color(this, i)) {
		SbColor attrColor;
		if (this->getPointColor(i, attrColor))
		    glColor3f(attrColor[0], attrColor[1], attrColor[2]);
		else
		    set_vlist_default_gl_color(this);
	    }
	    vlist_gl_vertex_at(this, geom, i);
	}
    }
    glEnd();
    glPopAttrib();
}

void
SoBRLVListShape::computeBBox(SoAction *UNUSED(action), SbBox3f &box, SbVec3f &center)
{
    box.makeEmpty();
    if (!this->visible.getValue()) {
	center = SbVec3f(0.0f, 0.0f, 0.0f);
	return;
    }

    const SoBRLVListShape *geom = this->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++) {
	box.extendBy(geom->point[i]);
	float scale = 0.0f;
	if (this->getPointScale(i, scale) && scale > 0.0f) {
	    const SbVec3f radius(scale, scale, scale);
	    box.extendBy(geom->point[i] - radius);
	    box.extendBy(geom->point[i] + radius);
	}
    }

    center = box.isEmpty() ? SbVec3f(0.0f, 0.0f, 0.0f) : box.getCenter();
}

void
SoBRLVListShape::generatePrimitives(SoAction *action)
{
    if (!this->visible.getValue() || !this->selectable.getValue())
	return;

    SoPrimitiveVertex v0;
    SoPrimitiveVertex v1;
    SoPointDetail pointDetail;
    SbVec3f last;
    SbBool haveLast = FALSE;
    int segmentIndex = 0;
    const SoBRLVListShape *geom = this->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    v0.setNormal(0.0f, 0.0f, 1.0f);
    v1.setNormal(0.0f, 0.0f, 1.0f);
    v0.setDetail(&pointDetail);
    v1.setDetail(&pointDetail);

    for (int i = 0; i < n; i++) {
	switch (geom->command[i]) {
	    case MOVE:
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast) {
		    pointDetail.setCoordinateIndex(segmentIndex);
		    v0.setPoint(last);
		    v1.setPoint(geom->point[i]);
		    this->invokeLineSegmentCallbacks(action, &v0, &v1);
		    segmentIndex++;
		}
		last = geom->point[i];
		haveLast = TRUE;
		break;
	    case POINT:
		pointDetail.setCoordinateIndex(i);
		v0.setPoint(geom->point[i]);
		{
		    SbVec3f normal;
		    if (this->getPointNormal(i, normal))
			v0.setNormal(normal);
		    else
			v0.setNormal(0.0f, 0.0f, 1.0f);
		}
		this->invokePointCallbacks(action, &v0);
		break;
	    default:
		break;
	}
    }
}

SoDetail *
SoBRLVListShape::createLineSegmentDetail(SoRayPickAction *UNUSED(action),
	const SoPrimitiveVertex *v1,
	const SoPrimitiveVertex *UNUSED(v2),
	SoPickedPoint *UNUSED(pp))
{
    SoBRLPickDetail *detail = new SoBRLPickDetail;
    detail->setPath(this->sourcePath.getValue());
    detail->setSourceInstanceKey(this->ownerSourceInstanceKey.getValue());
    detail->setSourceName(this->sourceName.getValue());
    detail->setSourceType(this->sourceType.getValue());
    detail->setSourceId(this->sourceId.getValue());
    detail->setRegionId(this->regionId.getValue());
    detail->setAirCode(this->airCode.getValue());
    detail->setMaterialId(this->materialId.getValue());
    detail->setLos(this->los.getValue());
    detail->setMaterialColor(this->materialColorValid.getValue(),
			     this->materialColor.getValue());
    detail->setMaterialShader(this->materialShader.getValue());
    detail->setPrimitive(SoBRLPickDetail::LINE_SEGMENT, -1);
    detail->setEditIntent(this->editIntentId.getValue(),
			  this->editIntentRole.getValue());

    const SoDetail *vertexDetail = v1 ? v1->getDetail() : NULL;
    if (vertexDetail && vertexDetail->isOfType(SoPointDetail::getClassTypeId())) {
	const SoPointDetail *pointDetail = static_cast<const SoPointDetail *>(vertexDetail);
	detail->setPrimitive(SoBRLPickDetail::LINE_SEGMENT, pointDetail->getCoordinateIndex());
    }
    if (v1)
	detail->setModelPoint(v1->getPoint());
    return detail;
}

SoDetail *
SoBRLVListShape::createPointDetail(SoRayPickAction *UNUSED(action),
				   const SoPrimitiveVertex *v,
				   SoPickedPoint *UNUSED(pp))
{
    SoBRLPickDetail *detail = new SoBRLPickDetail;
    detail->setPath(this->sourcePath.getValue());
    detail->setSourceInstanceKey(this->ownerSourceInstanceKey.getValue());
    detail->setSourceName(this->sourceName.getValue());
    detail->setSourceType(this->sourceType.getValue());
    detail->setSourceId(this->sourceId.getValue());
    detail->setRegionId(this->regionId.getValue());
    detail->setAirCode(this->airCode.getValue());
    detail->setMaterialId(this->materialId.getValue());
    detail->setLos(this->los.getValue());
    detail->setMaterialColor(this->materialColorValid.getValue(),
			     this->materialColor.getValue());
    detail->setMaterialShader(this->materialShader.getValue());
    detail->setPrimitive(SoBRLPickDetail::POINT, -1);
    detail->setEditIntent(this->editIntentId.getValue(),
			  this->editIntentRole.getValue());

    const SoDetail *vertexDetail = v ? v->getDetail() : NULL;
    if (vertexDetail && vertexDetail->isOfType(SoPointDetail::getClassTypeId())) {
	const SoPointDetail *pointDetail = static_cast<const SoPointDetail *>(vertexDetail);
	detail->setPrimitive(SoBRLPickDetail::POINT, pointDetail->getCoordinateIndex());
    }
    if (v)
	detail->setModelPoint(v->getPoint());
    return detail;
}
