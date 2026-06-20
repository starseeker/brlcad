/*                 M E A S U R E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/measure_action.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

#include <float.h>
#include <math.h>
#include <vector>

SO_ACTION_SOURCE(SoBRLMeasureAction);

static SbBool
measure_source_full_detail_result_valid(const BRLObolSourceMeshRequest &sourceRequest,
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

struct measure_segment_record {
    SbString path;
    int primitiveIndex;
    SbVec3f a;
    SbVec3f b;
};

static SbVec3f
closest_point_on_segment(const SbVec3f &query, const SbVec3f &a, const SbVec3f &b)
{
    SbVec3f ab = b - a;
    float len2 = ab.sqrLength();
    if (len2 <= 0.0f)
	return a;

    float t = (query - a).dot(ab) / len2;
    if (t < 0.0f)
	t = 0.0f;
    if (t > 1.0f)
	t = 1.0f;
    return a + ab * t;
}

static SbVec3f
closest_point_on_triangle(const SbVec3f &query,
	const SbVec3f &a,
	const SbVec3f &b,
	const SbVec3f &c)
{
    SbVec3f ab = b - a;
    SbVec3f ac = c - a;
    SbVec3f ap = query - a;
    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f)
	return a;

    SbVec3f bp = query - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3)
	return b;

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
	float v = d1 / (d1 - d3);
	return a + ab * v;
    }

    SbVec3f cp = query - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6)
	return c;

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
	float w = d2 / (d2 - d6);
	return a + ac * w;
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
	float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
	return b + (c - b) * w;
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

static float
clamp_float(float value, float minValue, float maxValue)
{
    if (value < minValue)
	return minValue;
    if (value > maxValue)
	return maxValue;
    return value;
}

static SbBool
same_point(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length() <= 1.0e-5f;
}

static SbBool
shared_segment_vertex(const measure_segment_record &sa,
	const measure_segment_record &sb,
	SbVec3f &shared,
	SbVec3f &otherA,
	SbVec3f &otherB)
{
    if (same_point(sa.a, sb.a)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.b;
	return TRUE;
    }
    if (same_point(sa.a, sb.b)) {
	shared = sa.a;
	otherA = sa.b;
	otherB = sb.a;
	return TRUE;
    }
    if (same_point(sa.b, sb.a)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.b;
	return TRUE;
    }
    if (same_point(sa.b, sb.b)) {
	shared = sa.b;
	otherA = sa.a;
	otherB = sb.a;
	return TRUE;
    }
    return FALSE;
}

static SbBool
segment_angle_degrees(const SbVec3f &shared,
	const SbVec3f &otherA,
	const SbVec3f &otherB,
	float &angleDegrees)
{
    SbVec3f va = otherA - shared;
    SbVec3f vb = otherB - shared;
    float lenA = va.length();
    float lenB = vb.length();
    if (lenA <= 0.0f || lenB <= 0.0f)
	return FALSE;

    float cosAngle = clamp_float(va.dot(vb) / (lenA * lenB), -1.0f, 1.0f);
    angleDegrees = acosf(cosAngle) * (180.0f / 3.14159265358979323846f);
    return TRUE;
}

SoBRLMeasureAction::SoBRLMeasureAction(void) :
    queryPoint(0.0f, 0.0f, 0.0f),
    nearestPoint(0.0f, 0.0f, 0.0f),
    anglePoint(0.0f, 0.0f, 0.0f),
    nearestPath(""),
    anglePath(""),
    totalLength(0.0f),
    surfaceArea(0.0f),
    nearestDistance(FLT_MAX),
    angleDegrees(0.0f),
    angleDistance(FLT_MAX),
    shapeCount(0),
    segmentCount(0),
    triangleCount(0),
    nearestPrimitiveIndex(-1),
    anglePrimitiveIndexA(-1),
    anglePrimitiveIndexB(-1),
    nearestPrimitiveKind(NONE),
    coordinateSpace(WORLD_SPACE),
    selectionFilter(ALL_SELECTION),
    highlightFilter(ALL_HIGHLIGHT),
    geometryPolicy(SoBRLMeasureAction::FULL_DETAIL),
    skippedLodDisplayMeshCount(0),
    haveQueryPoint(FALSE),
    haveNearestPrimitive(FALSE),
    haveAngle(FALSE)
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeasureAction);
    this->bounds.makeEmpty();
}

SoBRLMeasureAction::~SoBRLMeasureAction(void)
{
}

void
SoBRLMeasureAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLMeasureAction, SoAction);
    SO_ENABLE(SoBRLMeasureAction, SoModelMatrixElement);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLMeasureAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLVListShape, SoBRLMeasureAction::vlistShapeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape, SoBRLMeasureAction::meshShapeAction);
}

void
SoBRLMeasureAction::setQueryPoint(const SbVec3f &point)
{
    this->queryPoint = point;
    this->haveQueryPoint = TRUE;
}

void
SoBRLMeasureAction::clearQueryPoint(void)
{
    this->haveQueryPoint = FALSE;
}

void
SoBRLMeasureAction::setCoordinateSpace(CoordinateSpace space)
{
    this->coordinateSpace = space;
}

SoBRLMeasureAction::CoordinateSpace
SoBRLMeasureAction::getCoordinateSpace(void) const
{
    return this->coordinateSpace;
}

void
SoBRLMeasureAction::setSelectionFilter(SelectionFilter filter)
{
    this->selectionFilter = filter;
}

SoBRLMeasureAction::SelectionFilter
SoBRLMeasureAction::getSelectionFilter(void) const
{
    return this->selectionFilter;
}

void
SoBRLMeasureAction::setHighlightFilter(HighlightFilter filter)
{
    this->highlightFilter = filter;
}

SoBRLMeasureAction::HighlightFilter
SoBRLMeasureAction::getHighlightFilter(void) const
{
    return this->highlightFilter;
}

void
SoBRLMeasureAction::setGeometryPolicy(GeometryPolicy policy)
{
    this->geometryPolicy = (policy == SoBRLMeasureAction::DISPLAY_LEVEL) ?
	SoBRLMeasureAction::DISPLAY_LEVEL : SoBRLMeasureAction::FULL_DETAIL;
}

SoBRLMeasureAction::GeometryPolicy
SoBRLMeasureAction::getGeometryPolicy(void) const
{
    return this->geometryPolicy;
}

unsigned int
SoBRLMeasureAction::getSkippedLodDisplayMeshCount(void) const
{
    return this->skippedLodDisplayMeshCount;
}

int
SoBRLMeasureAction::getSourceBackedFullDetailRequestCount(void) const
{
    return static_cast<int>(this->sourceBackedFullDetailRequests.size());
}

const BRLObolSourceMeshRequest &
SoBRLMeasureAction::getSourceBackedFullDetailRequest(int index) const
{
    return this->sourceBackedFullDetailRequests.at(static_cast<size_t>(index));
}

SbBool
SoBRLMeasureAction::makeSourceBackedFullDetailLodRequest(int index,
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
SoBRLMeasureAction::consumeSourceBackedFullDetailResult(
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result)
{
    if (!measure_source_full_detail_result_valid(sourceRequest, result))
	return FALSE;
    if (!this->selectionAllows(sourceRequest.selected) ||
	    !this->highlightAllows(sourceRequest.highlighted))
	return TRUE;

    size_t faceCount = result.mesh.coordIndex.size() / 3;
    SbBool measuredShape = FALSE;
    for (size_t i = 0; i < faceCount; i++) {
	int ia = result.mesh.coordIndex[i * 3];
	int ib = result.mesh.coordIndex[i * 3 + 1];
	int ic = result.mesh.coordIndex[i * 3 + 2];
	if (ia < 0 || ib < 0 || ic < 0 ||
		static_cast<size_t>(ia) >= result.mesh.points.size() ||
		static_cast<size_t>(ib) >= result.mesh.points.size() ||
		static_cast<size_t>(ic) >= result.mesh.points.size())
	    return FALSE;

	SbVec3f pointA = this->pointForCoordinateSpace(
		sourceRequest.localToWorld,
		result.mesh.points[static_cast<size_t>(ia)]);
	SbVec3f pointB = this->pointForCoordinateSpace(
		sourceRequest.localToWorld,
		result.mesh.points[static_cast<size_t>(ib)]);
	SbVec3f pointC = this->pointForCoordinateSpace(
		sourceRequest.localToWorld,
		result.mesh.points[static_cast<size_t>(ic)]);
	measuredShape = TRUE;
	this->measureTriangle(sourceRequest.path, static_cast<int>(i),
		pointA, pointB, pointC);
    }

    if (measuredShape)
	this->shapeCount++;

    return TRUE;
}

int
SoBRLMeasureAction::submitSourceBackedFullDetailRequests(
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
SoBRLMeasureAction::consumeSourceBackedFullDetailResults(
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
	    if (this->consumeSourceBackedFullDetailResult(
		    this->sourceBackedFullDetailRequests[i], results[j])) {
		used[j] = TRUE;
		consumed++;
	    }
	    break;
	}
    }

    return consumed;
}

SbBool
SoBRLMeasureAction::hasSegments(void) const
{
    return this->segmentCount > 0;
}

int
SoBRLMeasureAction::getShapeCount(void) const
{
    return this->shapeCount;
}

int
SoBRLMeasureAction::getSegmentCount(void) const
{
    return this->segmentCount;
}

SbBool
SoBRLMeasureAction::hasFaces(void) const
{
    return this->triangleCount > 0;
}

int
SoBRLMeasureAction::getTriangleCount(void) const
{
    return this->triangleCount;
}

float
SoBRLMeasureAction::getSurfaceArea(void) const
{
    return this->surfaceArea;
}

float
SoBRLMeasureAction::getTotalLength(void) const
{
    return this->totalLength;
}

const SbBox3f &
SoBRLMeasureAction::getBounds(void) const
{
    return this->bounds;
}

SbBool
SoBRLMeasureAction::hasNearestSegment(void) const
{
    return this->haveNearestPrimitive && this->nearestPrimitiveKind == LINE_SEGMENT;
}

SbBool
SoBRLMeasureAction::hasNearestPrimitive(void) const
{
    return this->haveNearestPrimitive;
}

SoBRLMeasureAction::PrimitiveKind
SoBRLMeasureAction::getNearestPrimitiveKind(void) const
{
    return this->nearestPrimitiveKind;
}

const SbString &
SoBRLMeasureAction::getNearestPath(void) const
{
    return this->nearestPath;
}

int
SoBRLMeasureAction::getNearestPrimitiveIndex(void) const
{
    return this->nearestPrimitiveIndex;
}

const SbVec3f &
SoBRLMeasureAction::getNearestPoint(void) const
{
    return this->nearestPoint;
}

float
SoBRLMeasureAction::getNearestDistance(void) const
{
    return this->nearestDistance;
}

SbBool
SoBRLMeasureAction::hasAngle(void) const
{
    return this->haveAngle;
}

float
SoBRLMeasureAction::getAngleDegrees(void) const
{
    return this->angleDegrees;
}

const SbVec3f &
SoBRLMeasureAction::getAnglePoint(void) const
{
    return this->anglePoint;
}

const SbString &
SoBRLMeasureAction::getAnglePath(void) const
{
    return this->anglePath;
}

int
SoBRLMeasureAction::getAnglePrimitiveIndexA(void) const
{
    return this->anglePrimitiveIndexA;
}

int
SoBRLMeasureAction::getAnglePrimitiveIndexB(void) const
{
    return this->anglePrimitiveIndexB;
}

float
SoBRLMeasureAction::getAngleDistance(void) const
{
    return this->angleDistance;
}

void
SoBRLMeasureAction::beginTraversal(SoNode *node)
{
    this->resetResults();
    this->traverse(node);
}

void
SoBRLMeasureAction::nodeAction(SoAction *action, SoNode *node)
{
    node->doAction(action);
}

void
SoBRLMeasureAction::vlistShapeAction(SoAction *action, SoNode *node)
{
    SoBRLMeasureAction *measureAction = static_cast<SoBRLMeasureAction *>(action);
    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    if (!shape->visible.getValue())
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    std::vector<measure_segment_record> measuredSegments;

    int localSegmentCount = shape->getSegmentCount();
    SbBool measuredShape = FALSE;

    for (int i = 0; i < localSegmentCount; i++) {
	SbVec3f a;
	SbVec3f b;
	if (!shape->getSegment(i, a, b))
	    continue;
	if (!measureAction->selectionAllows(shape->isPrimitiveSelected(i)))
	    continue;
	if (!measureAction->highlightAllows(shape->isPrimitiveHighlighted(i)))
	    continue;

	SbVec3f pointA = measureAction->pointForCoordinateSpace(localToWorld, a);
	SbVec3f pointB = measureAction->pointForCoordinateSpace(localToWorld, b);

	measuredShape = TRUE;
	measureAction->measureSegment(shape->sourcePath.getValue(), i, pointA, pointB);
	measure_segment_record record;
	record.path = shape->sourcePath.getValue();
	record.primitiveIndex = i;
	record.a = pointA;
	record.b = pointB;
	measuredSegments.push_back(record);
    }

    for (size_t i = 0; i < measuredSegments.size(); i++) {
	for (size_t j = i + 1; j < measuredSegments.size(); j++) {
	    SbVec3f shared;
	    SbVec3f otherA;
	    SbVec3f otherB;
	    float degrees = 0.0f;
	    if (!shared_segment_vertex(measuredSegments[i], measuredSegments[j],
		    shared, otherA, otherB))
		continue;
	    if (!segment_angle_degrees(shared, otherA, otherB, degrees))
		continue;
	    measureAction->considerAngle(measuredSegments[i].path,
		    measuredSegments[i].primitiveIndex,
		    measuredSegments[j].primitiveIndex,
		    shared,
		    degrees);
	}
    }

    if (measuredShape)
	measureAction->shapeCount++;
}

void
SoBRLMeasureAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLMeasureAction *measureAction = static_cast<SoBRLMeasureAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    if (!shape->visible.getValue())
	return;
    const SbBool useFullDetail =
	(measureAction->geometryPolicy == SoBRLMeasureAction::FULL_DETAIL &&
	 shape->hasFullDetailMesh()) ? TRUE : FALSE;
    const SbBool useSourceBackedFullDetail =
	(measureAction->geometryPolicy == SoBRLMeasureAction::FULL_DETAIL &&
	 shape->needsSourceBackedFullDetail()) ? TRUE : FALSE;
    if (measureAction->geometryPolicy == SoBRLMeasureAction::FULL_DETAIL &&
	    shape->isLodDisplayActive())
	measureAction->skippedLodDisplayMeshCount++;
    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    if (useSourceBackedFullDetail) {
	measureAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
	return;
    }

    int localTriangleCount = useFullDetail ?
	shape->getFullDetailTriangleCount() : shape->getTriangleCount();
    SbBool measuredShape = FALSE;

    for (int i = 0; i < localTriangleCount; i++) {
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (useFullDetail) {
	    if (!shape->getFullDetailTriangle(i, a, b, c))
		continue;
	} else {
	    if (!shape->getTriangle(i, a, b, c))
		continue;
	}
	if (!measureAction->selectionAllows(shape->isPrimitiveSelected(i)))
	    continue;
	if (!measureAction->highlightAllows(shape->isPrimitiveHighlighted(i)))
	    continue;

	SbVec3f pointA = measureAction->pointForCoordinateSpace(localToWorld, a);
	SbVec3f pointB = measureAction->pointForCoordinateSpace(localToWorld, b);
	SbVec3f pointC = measureAction->pointForCoordinateSpace(localToWorld, c);

	measuredShape = TRUE;
	measureAction->measureTriangle(shape->sourcePath.getValue(), i, pointA, pointB, pointC);
    }

    if (measuredShape)
	measureAction->shapeCount++;
}

void
SoBRLMeasureAction::resetResults(void)
{
    this->nearestPoint = SbVec3f(0.0f, 0.0f, 0.0f);
    this->nearestPath = "";
    this->bounds.makeEmpty();
    this->totalLength = 0.0f;
    this->surfaceArea = 0.0f;
    this->nearestDistance = FLT_MAX;
    this->angleDegrees = 0.0f;
    this->angleDistance = FLT_MAX;
    this->shapeCount = 0;
    this->segmentCount = 0;
    this->triangleCount = 0;
    this->nearestPrimitiveIndex = -1;
    this->anglePrimitiveIndexA = -1;
    this->anglePrimitiveIndexB = -1;
    this->skippedLodDisplayMeshCount = 0;
    this->sourceBackedFullDetailRequests.clear();
    this->nearestPrimitiveKind = NONE;
    this->haveNearestPrimitive = FALSE;
    this->haveAngle = FALSE;
    this->anglePoint = SbVec3f(0.0f, 0.0f, 0.0f);
    this->anglePath = "";
}

void
SoBRLMeasureAction::appendSourceBackedFullDetailRequest(
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
SoBRLMeasureAction::measureSegment(const SbString &path, int primitiveIndex,
	const SbVec3f &a, const SbVec3f &b)
{
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
    this->totalLength += (b - a).length();
    this->segmentCount++;

    if (!this->haveQueryPoint)
	return;

    SbVec3f candidate = closest_point_on_segment(this->queryPoint, a, b);
    float dist = (this->queryPoint - candidate).length();
    if (dist >= this->nearestDistance)
	return;

    this->haveNearestPrimitive = TRUE;
    this->nearestPrimitiveKind = LINE_SEGMENT;
    this->nearestDistance = dist;
    this->nearestPoint = candidate;
    this->nearestPath = path;
    this->nearestPrimitiveIndex = primitiveIndex;
}

void
SoBRLMeasureAction::considerAngle(const SbString &path,
	int primitiveIndexA,
	int primitiveIndexB,
	const SbVec3f &point,
	float degrees)
{
    float dist = this->haveQueryPoint ? (this->queryPoint - point).length() : 0.0f;
    if (this->haveAngle) {
	if (!this->haveQueryPoint)
	    return;
	if (dist >= this->angleDistance)
	    return;
    }

    this->haveAngle = TRUE;
    this->anglePath = path;
    this->anglePoint = point;
    this->angleDegrees = degrees;
    this->angleDistance = dist;
    this->anglePrimitiveIndexA = primitiveIndexA;
    this->anglePrimitiveIndexB = primitiveIndexB;
}

SbVec3f
SoBRLMeasureAction::pointForCoordinateSpace(const SbMatrix &localToWorld,
	const SbVec3f &localPoint) const
{
    if (this->coordinateSpace == PATH_LOCAL_SPACE)
	return localPoint;

    SbVec3f ret;
    localToWorld.multVecMatrix(localPoint, ret);
    return ret;
}

void
SoBRLMeasureAction::measureTriangle(const SbString &path, int primitiveIndex,
	const SbVec3f &a, const SbVec3f &b, const SbVec3f &c)
{
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
    this->bounds.extendBy(c);
    this->surfaceArea += 0.5f * (b - a).cross(c - a).length();
    this->triangleCount++;

    if (!this->haveQueryPoint)
	return;

    SbVec3f candidate = closest_point_on_triangle(this->queryPoint, a, b, c);
    float dist = (this->queryPoint - candidate).length();
    if (dist >= this->nearestDistance)
	return;

    this->haveNearestPrimitive = TRUE;
    this->nearestPrimitiveKind = FACE;
    this->nearestDistance = dist;
    this->nearestPoint = candidate;
    this->nearestPath = path;
    this->nearestPrimitiveIndex = primitiveIndex;
}

SbBool
SoBRLMeasureAction::selectionAllows(SbBool selected) const
{
    if (this->selectionFilter == SELECTED_ONLY)
	return selected;
    if (this->selectionFilter == UNSELECTED_ONLY)
	return !selected;
    return TRUE;
}

SbBool
SoBRLMeasureAction::highlightAllows(SbBool highlighted) const
{
    if (this->highlightFilter == HIGHLIGHTED_ONLY)
	return highlighted;
    if (this->highlightFilter == UNHIGHLIGHTED_ONLY)
	return !highlighted;
    return TRUE;
}
