/*                 E X P O R T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/export_action.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

#include <vector>

SO_ACTION_SOURCE(SoBRLExportAction);

static SbBool
export_source_full_detail_result_valid(const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result)
{
    if (result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	    result.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	    !result.mesh.isValid())
	return FALSE;

    size_t faceCount = result.mesh.coordIndex.size() / 3;
    if ((sourceRequest.faceCount != 0 &&
		sourceRequest.faceCount != static_cast<uint64_t>(faceCount)) ||
	    (sourceRequest.pointCount != 0 &&
		sourceRequest.pointCount !=
		static_cast<uint64_t>(result.mesh.points.size())))
	return FALSE;

    return TRUE;
}

static int
export_source_mesh_face_index(const BRLObolLodMeshPayload &mesh,
	size_t faceSlot, size_t faceCount)
{
    if (mesh.faceIndex.size() == faceCount)
	return static_cast<int>(mesh.faceIndex[faceSlot]);
    return static_cast<int>(faceSlot);
}

static int
export_source_mesh_vertex_index(const BRLObolLodMeshPayload &mesh,
	int localIndex)
{
    if (localIndex < 0)
	return localIndex;
    if (mesh.vertexIndex.size() == mesh.points.size())
	return static_cast<int>(
		mesh.vertexIndex[static_cast<size_t>(localIndex)]);
    return localIndex;
}

static SbVec3f
export_transform_point_normal(const SbMatrix &localToWorld,
	const SbVec3f &localNormal)
{
    SbMatrix normalMatrix = localToWorld.inverse().transpose();
    SbVec3f worldNormal;
    normalMatrix.multDirMatrix(localNormal, worldNormal);
    if (worldNormal.length() > 0.0f)
	worldNormal.normalize();
    return worldNormal;
}

static float
export_transform_point_scale(const SbMatrix &localToWorld,
	const SbVec3f &localPoint,
	const SbVec3f &worldPoint,
	float localScale)
{
    if (localScale <= 0.0f)
	return localScale;

    const float offset[2] = {-localScale, localScale};
    float worldScale = 0.0f;

    for (int xi = 0; xi < 2; xi++) {
	for (int yi = 0; yi < 2; yi++) {
	    for (int zi = 0; zi < 2; zi++) {
		SbVec3f localCorner(
			localPoint[0] + offset[xi],
			localPoint[1] + offset[yi],
			localPoint[2] + offset[zi]);
		SbVec3f worldCorner;
		localToWorld.multVecMatrix(localCorner, worldCorner);
		const SbVec3f delta = worldCorner - worldPoint;
		for (int axis = 0; axis < 3; axis++) {
		    const float radius = delta[axis] < 0.0f ?
			-delta[axis] : delta[axis];
		    if (radius > worldScale)
			worldScale = radius;
		}
	    }
	}
    }

    return worldScale;
}

template <typename Record>
static void
export_reserve_records(std::vector<Record> &records, size_t additionalCount)
{
    if (additionalCount == 0)
	return;

    const size_t required = records.size() + additionalCount;
    if (required <= records.capacity())
	return;

    size_t grown = records.capacity() > 0 ?
	records.capacity() * 2 : additionalCount;
    if (grown < required)
	grown = required;
    records.reserve(grown);
}

SoBRLExportAction::SoBRLExportAction(void) :
    lineCount(0),
    pointCount(0),
    triangleCount(0),
    geometryPolicy(FULL_DETAIL),
    recordStorageEnabled(TRUE),
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
    return this->lineCount;
}

const SoBRLExportAction::LineRecord &
SoBRLExportAction::getLine(int index) const
{
    return this->lines.at(static_cast<size_t>(index));
}

int
SoBRLExportAction::getPointCount(void) const
{
    return this->pointCount;
}

const SoBRLExportAction::PointRecord &
SoBRLExportAction::getPoint(int index) const
{
    return this->points.at(static_cast<size_t>(index));
}

int
SoBRLExportAction::getTriangleCount(void) const
{
    return this->triangleCount;
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
SoBRLExportAction::setRecordStorageEnabled(SbBool enabled)
{
    this->recordStorageEnabled = enabled ? TRUE : FALSE;
}

SbBool
SoBRLExportAction::isRecordStorageEnabled(void) const
{
    return this->recordStorageEnabled;
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

int
SoBRLExportAction::getSourceBackedFullDetailRequestCount(void) const
{
    return static_cast<int>(this->sourceBackedFullDetailRequests.size());
}

const BRLObolSourceMeshRequest &
SoBRLExportAction::getSourceBackedFullDetailRequest(int index) const
{
    return this->sourceBackedFullDetailRequests.at(static_cast<size_t>(index));
}

SbBool
SoBRLExportAction::makeSourceBackedFullDetailLodRequest(int index,
	BRLObolLodRequest &request,
	const BRLObolLodRequest *templateRequest) const
{
    if (index < 0 ||
	    static_cast<size_t>(index) >=
	    this->sourceBackedFullDetailRequests.size())
	return FALSE;

    return brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	request, this->sourceBackedFullDetailRequests[static_cast<size_t>(index)],
	templateRequest);
}

SbBool
SoBRLExportAction::appendSourceBackedFullDetailResult(
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result)
{
    if (!export_source_full_detail_result_valid(sourceRequest, result))
	return FALSE;

    size_t faceCount = result.mesh.coordIndex.size() / 3;
    if (this->recordStorageEnabled)
	export_reserve_records(this->triangles, faceCount);
    for (size_t i = 0; i < faceCount; i++) {
	int ia = result.mesh.coordIndex[i * 3];
	int ib = result.mesh.coordIndex[i * 3 + 1];
	int ic = result.mesh.coordIndex[i * 3 + 2];
	if (ia < 0 || ib < 0 || ic < 0 ||
		static_cast<size_t>(ia) >= result.mesh.points.size() ||
		static_cast<size_t>(ib) >= result.mesh.points.size() ||
		static_cast<size_t>(ic) >= result.mesh.points.size())
	    return FALSE;

	SbVec3f worldA;
	SbVec3f worldB;
	SbVec3f worldC;
	sourceRequest.localToWorld.multVecMatrix(
		result.mesh.points[static_cast<size_t>(ia)], worldA);
	sourceRequest.localToWorld.multVecMatrix(
		result.mesh.points[static_cast<size_t>(ib)], worldB);
	sourceRequest.localToWorld.multVecMatrix(
		result.mesh.points[static_cast<size_t>(ic)], worldC);
	const int sourceIa = export_source_mesh_vertex_index(result.mesh, ia);
	const int sourceIb = export_source_mesh_vertex_index(result.mesh, ib);
	const int sourceIc = export_source_mesh_vertex_index(result.mesh, ic);

	if (!this->recordStorageEnabled) {
	    this->appendTriangleSummary(worldA, worldB, worldC);
	    continue;
	}

	this->appendTriangle(sourceRequest.path, sourceRequest.sourceName,
		sourceRequest.sourceType, sourceRequest.sourceId,
		sourceRequest.regionId, sourceRequest.airCode,
		sourceRequest.materialId, sourceRequest.los,
		sourceRequest.materialColorValid, sourceRequest.materialColor,
		sourceRequest.materialShader,
		export_source_mesh_face_index(result.mesh, i, faceCount),
		sourceIa, sourceIb, sourceIc, sourceRequest.selected,
		sourceRequest.highlighted,
		sourceRequest.ghosted, sourceRequest.hiddenLine,
		sourceRequest.editEmphasis, sourceRequest.editIntentId,
		sourceRequest.editIntentRole, sourceRequest.lodPolicy,
		sourceRequest.lodAvailable, sourceRequest.lodActiveLevel,
		result.counts.faceCount > 0 ?
		    static_cast<uint32_t>(result.counts.faceCount) :
		    static_cast<uint32_t>(faceCount),
		result.counts.pointCount > 0 ?
		    static_cast<uint32_t>(result.counts.pointCount) :
		    static_cast<uint32_t>(result.mesh.points.size()),
		static_cast<uint32_t>(result.counts.originalPointCount),
		static_cast<uint32_t>(result.counts.normalCount),
		result.hasSnappedPoints, result.hasNormals,
		result.bounds.isEmpty() ? sourceRequest.lodBoundsMin :
		    result.bounds.getMin(),
		result.bounds.isEmpty() ? sourceRequest.lodBoundsMax :
		    result.bounds.getMax(),
		sourceRequest.colorOverride, sourceRequest.color,
		worldA, worldB, worldC);
    }

    return TRUE;
}

int
SoBRLExportAction::submitSourceBackedFullDetailRequests(
	BRLObolLodService *service, uint64_t generation, struct db_i *dbip,
	const BRLObolLodRequest *templateRequest,
	uint64_t maxFullDetailFaceCount,
	uint64_t maxFullDetailPointCount) const
{
    int submitted = 0;
    for (size_t i = 0; i < this->sourceBackedFullDetailRequests.size(); i++) {
	if (brlobol_lod_submit_rt_source_full_detail_request(service,
		generation, this->sourceBackedFullDetailRequests[i], dbip,
		templateRequest, maxFullDetailFaceCount,
		maxFullDetailPointCount) != 0)
	    submitted++;
    }
    return submitted;
}

int
SoBRLExportAction::consumeSourceBackedFullDetailResults(
	const std::vector<BRLObolLodResult> &results,
	const BRLObolLodRequest *templateRequest)
{
    std::vector<SbBool> used(results.size(), FALSE);
    int consumed = 0;

    for (size_t i = 0; i < this->sourceBackedFullDetailRequests.size(); i++) {
	BRLObolLodRequest expected;
	if (!brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
		expected, this->sourceBackedFullDetailRequests[i],
		templateRequest))
	    continue;

	for (size_t j = 0; j < results.size(); j++) {
	    if (used[j] ||
		    !brlobol_lod_result_matches_request(results[j], expected))
		continue;
	    if (this->appendSourceBackedFullDetailResult(
		    this->sourceBackedFullDetailRequests[i], results[j])) {
		used[j] = TRUE;
		consumed++;
	    }
	    break;
	}
    }

    return consumed;
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
    if (exportAction->recordStorageEnabled) {
	export_reserve_records(exportAction->lines,
		static_cast<size_t>(shape->getSegmentCount()));
	export_reserve_records(exportAction->points,
		static_cast<size_t>(shape->getPointPrimitiveCount()));
    }

    SbVec3f last;
    SbBool haveLast = FALSE;
    int segmentIndex = 0;
    int n = shape->point.getNum();
    if (shape->command.getNum() < n)
	n = shape->command.getNum();

    for (int i = 0; i < n; i++) {
	if (shape->command[i] == SoBRLVListShape::MOVE) {
	    last = shape->point[i];
	    haveLast = TRUE;
	    continue;
	}

	if (shape->command[i] == SoBRLVListShape::DRAW) {
	    if (haveLast) {
		SbVec3f worldA;
		SbVec3f worldB;
		localToWorld.multVecMatrix(last, worldA);
		localToWorld.multVecMatrix(shape->point[i], worldB);

		if (!exportAction->recordStorageEnabled) {
		    exportAction->appendLineSummary(worldA, worldB);
		} else {
		    exportAction->appendLine(shape->sourcePath.getValue(),
			    shape->sourceName.getValue(), shape->sourceType.getValue(),
			    shape->sourceId.getValue(),
			    shape->regionId.getValue(), shape->airCode.getValue(),
			    shape->materialId.getValue(), shape->los.getValue(),
			    shape->materialColorValid.getValue(),
			    shape->materialColor.getValue(),
			    shape->materialShader.getValue(), segmentIndex,
			    shape->isPrimitiveSelected(segmentIndex),
			    shape->isPrimitiveHighlighted(segmentIndex),
			    shape->ghosted.getValue(), shape->hiddenLine.getValue(),
			    shape->editEmphasis.getValue(),
			    shape->editIntentId.getValue(),
			    shape->editIntentRole.getValue(),
			    shape->lodPolicy.getValue(),
			    shape->colorOverride.getValue(),
			    shape->color.getValue(), worldA, worldB);
		}
		segmentIndex++;
	    }
	    last = shape->point[i];
	    haveLast = TRUE;
	    continue;
	}

	if (shape->command[i] == SoBRLVListShape::POINT) {
	    SbVec3f worldPoint;
	    localToWorld.multVecMatrix(shape->point[i], worldPoint);

	    float pointScale = 0.0f;
	    const int pointScaleValid =
		shape->getPointScale(i, pointScale) ? 1 : 0;
	    if (pointScaleValid)
		pointScale = export_transform_point_scale(localToWorld,
			shape->point[i], worldPoint, pointScale);

	    if (!exportAction->recordStorageEnabled) {
		exportAction->appendPointSummary(pointScaleValid, pointScale,
			worldPoint);
		continue;
	    }

	    SbColor pointColor(1.0f, 1.0f, 1.0f);
	    SbVec3f pointNormal(0.0f, 0.0f, 1.0f);
	    const int pointColorValid =
		shape->getPointColor(i, pointColor) ? 1 : 0;
	    const int pointNormalValid =
		shape->getPointNormal(i, pointNormal) ? 1 : 0;
	    if (pointNormalValid)
		pointNormal = export_transform_point_normal(localToWorld,
			pointNormal);

	    exportAction->appendPoint(shape->sourcePath.getValue(),
		    shape->sourceName.getValue(), shape->sourceType.getValue(),
		    shape->sourceId.getValue(),
		    shape->regionId.getValue(), shape->airCode.getValue(),
		    shape->materialId.getValue(), shape->los.getValue(),
		    shape->materialColorValid.getValue(),
		    shape->materialColor.getValue(),
		    shape->materialShader.getValue(), i,
		    shape->isPrimitiveSelected(i),
		    shape->isPrimitiveHighlighted(i),
		    shape->ghosted.getValue(), shape->hiddenLine.getValue(),
		    shape->editEmphasis.getValue(),
		    shape->editIntentId.getValue(),
		    shape->editIntentRole.getValue(),
		    shape->lodPolicy.getValue(),
		    shape->colorOverride.getValue(),
		    shape->color.getValue(),
		    pointColorValid, pointColor,
		    pointScaleValid, pointScale,
		    pointNormalValid, pointNormal,
		    worldPoint);
	}
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
	 shape->hasFullDetailMesh()) ? TRUE : FALSE;
    const SbBool useSourceBackedFullDetail =
	(exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	 shape->needsSourceBackedFullDetail()) ? TRUE : FALSE;
    if (exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	    shape->isLodDisplayActive())
	exportAction->skippedLodDisplayMeshCount++;
    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    if (useSourceBackedFullDetail) {
	exportAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
	return;
    }

    int triangleCount = useFullDetail ?
	shape->getFullDetailTriangleCount() : shape->getTriangleCount();
    if (exportAction->recordStorageEnabled)
	export_reserve_records(exportAction->triangles,
		static_cast<size_t>(triangleCount));
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

	if (!exportAction->recordStorageEnabled) {
	    exportAction->appendTriangleSummary(worldA, worldB, worldC);
	    continue;
	}

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
		shape->editEmphasis.getValue(),
		shape->editIntentId.getValue(),
		shape->editIntentRole.getValue(),
		shape->lodPolicy.getValue(),
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
    this->sourceBackedFullDetailRequests.clear();
    this->lineCount = 0;
    this->pointCount = 0;
    this->triangleCount = 0;
    this->bounds.makeEmpty();
    this->skippedLodDisplayMeshCount = 0;
}

void
SoBRLExportAction::appendSourceBackedFullDetailRequest(
	const SoBRLMeshShape *shape, const SbMatrix &localToWorld)
{
    if (!shape)
	return;

    BRLObolSourceMeshRequest request;
    if (shape->makeSourceMeshRequest(request)) {
	request.localToWorld = localToWorld;
	this->sourceBackedFullDetailRequests.push_back(request);
    }
}

void
SoBRLExportAction::appendLineSummary(const SbVec3f &a, const SbVec3f &b)
{
    this->lineCount++;
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
}

void
SoBRLExportAction::appendPointSummary(int pointScaleValid, float pointScale,
	const SbVec3f &point)
{
    this->pointCount++;
    this->bounds.extendBy(point);
    if (pointScaleValid && pointScale > 0.0f) {
	const SbVec3f radius(pointScale, pointScale, pointScale);
	this->bounds.extendBy(point - radius);
	this->bounds.extendBy(point + radius);
    }
}

void
SoBRLExportAction::appendTriangleSummary(const SbVec3f &a, const SbVec3f &b,
	const SbVec3f &c)
{
    this->triangleCount++;
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
    this->bounds.extendBy(c);
}

void
SoBRLExportAction::appendLine(const SbString &path, const SbString &sourceName,
	const SbString &sourceType, uint32_t sourceId,
	int regionId, int airCode, int materialId, int los,
	int materialColorValid, const SbColor &materialColor,
	const SbString &materialShader, int primitiveIndex,
	int selected, int highlighted, int ghosted,
	int hiddenLine, int editEmphasis,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	uint32_t lodPolicy,
	int colorOverride, const SbColor &color,
	const SbVec3f &a, const SbVec3f &b)
{
    if (!this->recordStorageEnabled) {
	this->appendLineSummary(a, b);
	return;
    }

    this->lineCount++;
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
    record.editIntentId = editIntentId;
    record.editIntentRole = editIntentRole;
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
	int hiddenLine, int editEmphasis,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	uint32_t lodPolicy,
	int colorOverride, const SbColor &color,
	int pointColorValid, const SbColor &pointColor,
	int pointScaleValid, float pointScale,
	int pointNormalValid, const SbVec3f &pointNormal,
	const SbVec3f &point)
{
    if (!this->recordStorageEnabled) {
	this->appendPointSummary(pointScaleValid, pointScale, point);
	return;
    }

    this->pointCount++;
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
    record.editIntentId = editIntentId;
    record.editIntentRole = editIntentRole;
    record.lodPolicy = lodPolicy;
    record.colorOverride = colorOverride ? 1 : 0;
    record.color = color;
    record.pointColorValid = pointColorValid ? 1 : 0;
    record.pointColor = pointColor;
    record.pointScaleValid = pointScaleValid ? 1 : 0;
    record.pointScale = pointScale;
    record.pointNormalValid = pointNormalValid ? 1 : 0;
    record.pointNormal = pointNormal;
    record.point = point;
    this->points.push_back(record);
    this->bounds.extendBy(point);
    if (record.pointScaleValid && record.pointScale > 0.0f) {
	const SbVec3f radius(record.pointScale, record.pointScale,
		record.pointScale);
	this->bounds.extendBy(point - radius);
	this->bounds.extendBy(point + radius);
    }
}

void
SoBRLExportAction::appendTriangle(const SbString &path,
	const SbString &sourceName, const SbString &sourceType,
	uint32_t sourceId, int regionId, int airCode, int materialId, int los,
	const int materialColorValid, const SbColor &materialColor,
	const SbString &materialShader, int primitiveIndex,
	int vertexIndexA, int vertexIndexB, int vertexIndexC,
	int selected, int highlighted, int ghosted,
	int hiddenLine, int editEmphasis,
	const SbString &editIntentId,
	const SbString &editIntentRole,
	uint32_t lodPolicy,
	int lodAvailable, int lodActiveLevel, uint32_t lodFaceCount,
	uint32_t lodPointCount, uint32_t lodOriginalPointCount,
	uint32_t lodNormalCount, int lodHasSnappedPoints,
	int lodHasNormals, const SbVec3f &lodBoundsMin,
	const SbVec3f &lodBoundsMax,
	const int colorOverride, const SbColor &color,
	const SbVec3f &a, const SbVec3f &b, const SbVec3f &c)
{
    if (!this->recordStorageEnabled) {
	this->appendTriangleSummary(a, b, c);
	return;
    }

    this->triangleCount++;
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
    record.editIntentId = editIntentId;
    record.editIntentRole = editIntentRole;
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
