/*                 E X P O R T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/export_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

SO_ACTION_SOURCE(SoBRLExportAction);

SoBRLExportAction::SoBRLExportAction(void) :
    geometryPolicy(FULL_DETAIL),
    skippedLodDisplayMeshCount(0)
{
    SO_ACTION_CONSTRUCTOR(SoBRLExportAction);
    this->bounds.makeEmpty();
}

SoBRLExportAction::~SoBRLExportAction(void)
{
}

void
SoBRLExportAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLExportAction, SoAction);
    SO_ENABLE(SoBRLExportAction, SoModelMatrixElement);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLExportAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLVListShape, SoBRLExportAction::vlistShapeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape, SoBRLExportAction::meshShapeAction);
}

int
SoBRLExportAction::getLineCount(void) const
{
    return static_cast<int>(this->lines.size());
}

const SoBRLExportAction::LineRecord &
SoBRLExportAction::getLine(int index) const
{
    return this->lines.at(static_cast<size_t>(index));
}

int
SoBRLExportAction::getPointCount(void) const
{
    return static_cast<int>(this->points.size());
}

const SoBRLExportAction::PointRecord &
SoBRLExportAction::getPoint(int index) const
{
    return this->points.at(static_cast<size_t>(index));
}

int
SoBRLExportAction::getTriangleCount(void) const
{
    return static_cast<int>(this->triangles.size());
}

const SoBRLExportAction::TriangleRecord &
SoBRLExportAction::getTriangle(int index) const
{
    return this->triangles.at(static_cast<size_t>(index));
}

const SbBox3f &
SoBRLExportAction::getBounds(void) const
{
    return this->bounds;
}

void
SoBRLExportAction::setGeometryPolicy(GeometryPolicy policy)
{
    this->geometryPolicy = (policy == SoBRLExportAction::DISPLAY_LEVEL) ?
	SoBRLExportAction::DISPLAY_LEVEL : SoBRLExportAction::FULL_DETAIL;
}

SoBRLExportAction::GeometryPolicy
SoBRLExportAction::getGeometryPolicy(void) const
{
    return this->geometryPolicy;
}

unsigned int
SoBRLExportAction::getSkippedLodDisplayMeshCount(void) const
{
    return this->skippedLodDisplayMeshCount;
}

void
SoBRLExportAction::beginTraversal(SoNode *node)
{
    this->resetResults();
    this->traverse(node);
}

void
SoBRLExportAction::nodeAction(SoAction *action, SoNode *node)
{
    node->doAction(action);
}

void
SoBRLExportAction::vlistShapeAction(SoAction *action, SoNode *node)
{
    SoBRLExportAction *exportAction = static_cast<SoBRLExportAction *>(action);
    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    if (!shape->visible.getValue())
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());

    for (int i = 0; i < shape->getSegmentCount(); i++) {
	SbVec3f a;
	SbVec3f b;
	if (!shape->getSegment(i, a, b))
	    continue;

	SbVec3f worldA;
	SbVec3f worldB;
	localToWorld.multVecMatrix(a, worldA);
	localToWorld.multVecMatrix(b, worldB);

	exportAction->appendLine(shape->sourcePath.getValue(),
		shape->sourceName.getValue(), shape->sourceType.getValue(),
		shape->sourceId.getValue(),
		shape->regionId.getValue(), shape->airCode.getValue(),
		shape->materialId.getValue(), shape->los.getValue(),
		shape->materialColorValid.getValue(),
		shape->materialColor.getValue(),
		shape->materialShader.getValue(), i,
		shape->isPrimitiveSelected(i), shape->isPrimitiveHighlighted(i),
		shape->ghosted.getValue(), shape->hiddenLine.getValue(),
		shape->editEmphasis.getValue(), shape->lodPolicy.getValue(),
		shape->colorOverride.getValue(),
		shape->color.getValue(), worldA, worldB);
    }

    for (int i = 0; i < shape->getPointPrimitiveCount(); i++) {
	SbVec3f point;
	int primitiveIndex = -1;
	if (!shape->getPointPrimitive(i, primitiveIndex, point))
	    continue;

	SbVec3f worldPoint;
	localToWorld.multVecMatrix(point, worldPoint);

	exportAction->appendPoint(shape->sourcePath.getValue(),
		shape->sourceName.getValue(), shape->sourceType.getValue(),
		shape->sourceId.getValue(),
		shape->regionId.getValue(), shape->airCode.getValue(),
		shape->materialId.getValue(), shape->los.getValue(),
		shape->materialColorValid.getValue(),
		shape->materialColor.getValue(),
		shape->materialShader.getValue(), primitiveIndex,
		shape->isPrimitiveSelected(primitiveIndex),
		shape->isPrimitiveHighlighted(primitiveIndex),
		shape->ghosted.getValue(), shape->hiddenLine.getValue(),
		shape->editEmphasis.getValue(), shape->lodPolicy.getValue(),
		shape->colorOverride.getValue(),
		shape->color.getValue(), worldPoint);
    }
}

void
SoBRLExportAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLExportAction *exportAction = static_cast<SoBRLExportAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    if (!shape->visible.getValue())
	return;

    const SbBool useFullDetail =
	(exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	 shape->isLodDisplayActive() && shape->hasFullDetailMesh()) ?
	TRUE : FALSE;
    if (exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	    shape->isLodDisplayActive())
	exportAction->skippedLodDisplayMeshCount++;
    if (exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	    shape->isLodDisplayActive() && !useFullDetail)
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());

    int triangleCount = useFullDetail ?
	shape->getFullDetailTriangleCount() : shape->getTriangleCount();
    for (int i = 0; i < triangleCount; i++) {
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	int ia = -1;
	int ib = -1;
	int ic = -1;
	if (useFullDetail) {
	    if (!shape->getFullDetailTriangle(i, a, b, c))
		continue;
	    if (!shape->getFullDetailTriangleVertexIndices(i, ia, ib, ic))
		continue;
	} else {
	    if (!shape->getTriangle(i, a, b, c))
		continue;
	    if (!shape->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
	}

	SbVec3f worldA;
	SbVec3f worldB;
	SbVec3f worldC;
	localToWorld.multVecMatrix(a, worldA);
	localToWorld.multVecMatrix(b, worldB);
	localToWorld.multVecMatrix(c, worldC);

	exportAction->appendTriangle(shape->sourcePath.getValue(),
		shape->sourceName.getValue(), shape->sourceType.getValue(),
		shape->sourceId.getValue(),
		shape->regionId.getValue(), shape->airCode.getValue(),
		shape->materialId.getValue(), shape->los.getValue(),
		shape->materialColorValid.getValue(),
		shape->materialColor.getValue(),
		shape->materialShader.getValue(), i,
		ia, ib, ic,
		shape->isPrimitiveSelected(i), shape->isPrimitiveHighlighted(i),
		shape->ghosted.getValue(), shape->hiddenLine.getValue(),
		shape->editEmphasis.getValue(), shape->lodPolicy.getValue(),
		shape->lodAvailable.getValue(),
		shape->lodActiveLevel.getValue(),
		shape->lodFaceCount.getValue(),
		shape->lodPointCount.getValue(),
		shape->lodOriginalPointCount.getValue(),
		shape->lodNormalCount.getValue(),
		shape->lodHasSnappedPoints.getValue(),
		shape->lodHasNormals.getValue(),
		shape->lodBoundsMin.getValue(),
		shape->lodBoundsMax.getValue(),
		shape->colorOverride.getValue(),
		shape->color.getValue(), worldA, worldB, worldC);
    }
}

void
SoBRLExportAction::resetResults(void)
{
    this->lines.clear();
    this->points.clear();
    this->triangles.clear();
    this->bounds.makeEmpty();
    this->skippedLodDisplayMeshCount = 0;
}

void
SoBRLExportAction::appendLine(const SbString &path, const SbString &sourceName,
	const SbString &sourceType, uint32_t sourceId,
	int regionId, int airCode, int materialId, int los,
	int materialColorValid, const SbColor &materialColor,
	const SbString &materialShader, int primitiveIndex,
	int selected, int highlighted, int ghosted,
	int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	int colorOverride, const SbColor &color,
	const SbVec3f &a, const SbVec3f &b)
{
    LineRecord record;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    record.regionId = regionId;
    record.airCode = airCode;
    record.materialId = materialId;
    record.los = los;
    record.materialColorValid = materialColorValid ? 1 : 0;
    record.materialColor = materialColor;
    record.materialShader = materialShader;
    record.primitiveIndex = primitiveIndex;
    record.selected = selected ? 1 : 0;
    record.highlighted = highlighted ? 1 : 0;
    record.ghosted = ghosted ? 1 : 0;
    record.hiddenLine = hiddenLine ? 1 : 0;
    record.editEmphasis = editEmphasis ? 1 : 0;
    record.lodPolicy = lodPolicy;
    record.colorOverride = colorOverride ? 1 : 0;
    record.color = color;
    record.a = a;
    record.b = b;
    this->lines.push_back(record);
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
}

void
SoBRLExportAction::appendPoint(const SbString &path, const SbString &sourceName,
	const SbString &sourceType, uint32_t sourceId,
	int regionId, int airCode, int materialId, int los,
	int materialColorValid, const SbColor &materialColor,
	const SbString &materialShader, int primitiveIndex,
	int selected, int highlighted, int ghosted,
	int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	int colorOverride, const SbColor &color,
	const SbVec3f &point)
{
    PointRecord record;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    record.regionId = regionId;
    record.airCode = airCode;
    record.materialId = materialId;
    record.los = los;
    record.materialColorValid = materialColorValid ? 1 : 0;
    record.materialColor = materialColor;
    record.materialShader = materialShader;
    record.primitiveIndex = primitiveIndex;
    record.selected = selected ? 1 : 0;
    record.highlighted = highlighted ? 1 : 0;
    record.ghosted = ghosted ? 1 : 0;
    record.hiddenLine = hiddenLine ? 1 : 0;
    record.editEmphasis = editEmphasis ? 1 : 0;
    record.lodPolicy = lodPolicy;
    record.colorOverride = colorOverride ? 1 : 0;
    record.color = color;
    record.point = point;
    this->points.push_back(record);
    this->bounds.extendBy(point);
}

void
SoBRLExportAction::appendTriangle(const SbString &path,
	const SbString &sourceName, const SbString &sourceType,
	uint32_t sourceId, int regionId, int airCode, int materialId, int los,
	const int materialColorValid, const SbColor &materialColor,
	const SbString &materialShader, int primitiveIndex,
	int vertexIndexA, int vertexIndexB, int vertexIndexC,
	int selected, int highlighted, int ghosted,
	int hiddenLine, int editEmphasis, uint32_t lodPolicy,
	int lodAvailable, int lodActiveLevel, uint32_t lodFaceCount,
	uint32_t lodPointCount, uint32_t lodOriginalPointCount,
	uint32_t lodNormalCount, int lodHasSnappedPoints,
	int lodHasNormals, const SbVec3f &lodBoundsMin,
	const SbVec3f &lodBoundsMax,
	const int colorOverride, const SbColor &color,
	const SbVec3f &a, const SbVec3f &b, const SbVec3f &c)
{
    TriangleRecord record;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    record.regionId = regionId;
    record.airCode = airCode;
    record.materialId = materialId;
    record.los = los;
    record.materialColorValid = materialColorValid ? 1 : 0;
    record.materialColor = materialColor;
    record.materialShader = materialShader;
    record.primitiveIndex = primitiveIndex;
    record.vertexIndexA = vertexIndexA;
    record.vertexIndexB = vertexIndexB;
    record.vertexIndexC = vertexIndexC;
    record.selected = selected ? 1 : 0;
    record.highlighted = highlighted ? 1 : 0;
    record.ghosted = ghosted ? 1 : 0;
    record.hiddenLine = hiddenLine ? 1 : 0;
    record.editEmphasis = editEmphasis ? 1 : 0;
    record.lodPolicy = lodPolicy;
    record.lodAvailable = lodAvailable ? 1 : 0;
    record.lodActiveLevel = lodActiveLevel;
    record.lodFaceCount = lodFaceCount;
    record.lodPointCount = lodPointCount;
    record.lodOriginalPointCount = lodOriginalPointCount;
    record.lodNormalCount = lodNormalCount;
    record.lodHasSnappedPoints = lodHasSnappedPoints ? 1 : 0;
    record.lodHasNormals = lodHasNormals ? 1 : 0;
    record.lodBoundsMin = lodBoundsMin;
    record.lodBoundsMax = lodBoundsMax;
    record.colorOverride = colorOverride ? 1 : 0;
    record.color = color;
    record.a = a;
    record.b = b;
    record.c = c;
    this->triangles.push_back(record);
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
    this->bounds.extendBy(c);
}
