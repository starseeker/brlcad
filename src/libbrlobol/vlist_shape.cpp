/*                   V L I S T _ S H A P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/pick_detail.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoPointDetail.h>
#include <Inventor/gl.h>

SO_NODE_SOURCE(SoBRLVListShape);

SoBRLVListShape::SoBRLVListShape(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLVListShape);

    SO_NODE_DEFINE_ENUM_VALUE(Command, MOVE);
    SO_NODE_DEFINE_ENUM_VALUE(Command, DRAW);
    SO_NODE_DEFINE_ENUM_VALUE(Command, POINT);

    SO_NODE_ADD_EMPTY_MFIELD(point);
    SO_NODE_ADD_EMPTY_MFIELD(command);
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
    SO_NODE_ADD_FIELD(regionId, (0));
    SO_NODE_ADD_FIELD(airCode, (0));
    SO_NODE_ADD_FIELD(materialId, (0));
    SO_NODE_ADD_FIELD(los, (0));
    SO_NODE_ADD_FIELD(materialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(materialColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialShader, (""));
    SO_NODE_ADD_FIELD(colorOverride, (FALSE));
    SO_NODE_ADD_FIELD(color, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(selectedColor, (SbColor(0.0f, 0.75f, 1.0f)));
    SO_NODE_ADD_FIELD(highlightedColor, (SbColor(1.0f, 1.0f, 0.0f)));
    SO_NODE_ADD_FIELD(ghostedColor, (SbColor(0.55f, 0.55f, 0.55f)));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selectable, (TRUE));
    SO_NODE_ADD_FIELD(selected, (FALSE));
    SO_NODE_ADD_FIELD(highlighted, (FALSE));
    SO_NODE_ADD_FIELD(ghosted, (FALSE));
    SO_NODE_ADD_FIELD(hiddenLine, (FALSE));
    SO_NODE_ADD_FIELD(editEmphasis, (FALSE));
    SO_NODE_ADD_FIELD(editIntentId, (""));
    SO_NODE_ADD_FIELD(editIntentRole, (""));
    SO_NODE_ADD_FIELD(lodPolicy, (0));
    SO_NODE_ADD_EMPTY_MFIELD(selectedPrimitive);
    SO_NODE_ADD_EMPTY_MFIELD(highlightedPrimitive);
}

SoBRLVListShape::~SoBRLVListShape(void)
{
}

void
SoBRLVListShape::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLVListShape, SoShape, "Shape");
}

void
SoBRLVListShape::setLineSet(const SbVec3f *points, const int32_t *commands, int count)
{
    this->point.setNum(0);
    this->command.setNum(0);
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
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    for (int i = 0; i < n; i++) {
	switch (this->command[i]) {
	    case MOVE:
		last = this->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast)
		    ret++;
		last = this->point[i];
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
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    for (int i = 0; i < n; i++) {
	switch (this->command[i]) {
	    case MOVE:
		last = this->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast) {
		    if (currentSegment == segmentIndex) {
			a = last;
			b = this->point[i];
			return TRUE;
		    }
		    currentSegment++;
		}
		last = this->point[i];
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
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    for (int i = 0; i < n; i++) {
	if (this->command[i] == POINT)
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
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    for (int i = 0; i < n; i++) {
	if (this->command[i] != POINT)
	    continue;
	if (currentPoint == pointIndex) {
	    primitiveIndex = i;
	    pointOut = this->point[i];
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
    if (primitiveIndex < 0 ||
	    primitiveIndex >= this->pointColorValid.getNum() ||
	    primitiveIndex >= this->pointColor.getNum() ||
	    !this->pointColorValid[primitiveIndex])
	return FALSE;

    colorOut = this->pointColor[primitiveIndex];
    return TRUE;
}

SbBool
SoBRLVListShape::getPointScale(int primitiveIndex, float &scaleOut) const
{
    scaleOut = 0.0f;
    if (primitiveIndex < 0 ||
	    primitiveIndex >= this->pointScaleValid.getNum() ||
	    primitiveIndex >= this->pointScale.getNum() ||
	    !this->pointScaleValid[primitiveIndex])
	return FALSE;

    scaleOut = this->pointScale[primitiveIndex];
    return TRUE;
}

SbBool
SoBRLVListShape::getPointNormal(int primitiveIndex, SbVec3f &normalOut) const
{
    normalOut = SbVec3f(0.0f, 0.0f, 1.0f);
    if (primitiveIndex < 0 ||
	    primitiveIndex >= this->pointNormalValid.getNum() ||
	    primitiveIndex >= this->pointNormal.getNum() ||
	    !this->pointNormalValid[primitiveIndex])
	return FALSE;

    normalOut = this->pointNormal[primitiveIndex];
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
set_vlist_gl_color(SoBRLVListShape *shape, int primitiveIndex)
{
    (void)set_vlist_gl_state_color(shape, primitiveIndex);
}

void
SoBRLVListShape::GLRender(SoGLRenderAction *action)
{
    if (!this->visible.getValue() || !this->shouldGLRender(action))
	return;

    glPushAttrib(GL_CURRENT_BIT);
    set_vlist_gl_color(this, -1);

    SbVec3f last;
    SbBool haveLast = FALSE;
    int segmentIndex = 0;
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    glBegin(GL_LINES);
    for (int i = 0; i < n; i++) {
	switch (this->command[i]) {
	    case MOVE:
		last = this->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast) {
		    set_vlist_gl_color(this, segmentIndex);
		    glVertex3f(last[0], last[1], last[2]);
		    glVertex3f(this->point[i][0], this->point[i][1], this->point[i][2]);
		    segmentIndex++;
		}
		last = this->point[i];
		haveLast = TRUE;
		break;
	    default:
		break;
	}
    }
    glEnd();

    glBegin(GL_POINTS);
    for (int i = 0; i < n; i++) {
	if (this->command[i] == POINT) {
	    if (!set_vlist_gl_state_color(this, i)) {
		SbColor attrColor;
		if (this->getPointColor(i, attrColor))
		    glColor3f(attrColor[0], attrColor[1], attrColor[2]);
	    }
	    glVertex3f(this->point[i][0], this->point[i][1], this->point[i][2]);
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

    for (int i = 0; i < this->point.getNum(); i++) {
	box.extendBy(this->point[i]);
	float scale = 0.0f;
	if (this->getPointScale(i, scale) && scale > 0.0f) {
	    const SbVec3f radius(scale, scale, scale);
	    box.extendBy(this->point[i] - radius);
	    box.extendBy(this->point[i] + radius);
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
    int n = this->point.getNum();
    if (this->command.getNum() < n)
	n = this->command.getNum();

    v0.setNormal(0.0f, 0.0f, 1.0f);
    v1.setNormal(0.0f, 0.0f, 1.0f);
    v0.setDetail(&pointDetail);
    v1.setDetail(&pointDetail);

    for (int i = 0; i < n; i++) {
	switch (this->command[i]) {
	    case MOVE:
		last = this->point[i];
		haveLast = TRUE;
		break;
	    case DRAW:
		if (haveLast) {
		    pointDetail.setCoordinateIndex(segmentIndex);
		    v0.setPoint(last);
		    v1.setPoint(this->point[i]);
		    this->invokeLineSegmentCallbacks(action, &v0, &v1);
		    segmentIndex++;
		}
		last = this->point[i];
		haveLast = TRUE;
		break;
	    case POINT:
		pointDetail.setCoordinateIndex(i);
		v0.setPoint(this->point[i]);
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
