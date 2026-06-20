/*                    S N A P _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/snap_action.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SbBox.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

#include <float.h>
#include <math.h>
#include <vector>

SO_ACTION_SOURCE(SoBRLSnapAction);

static SbBool
snap_source_full_detail_result_valid(const BRLObolSourceMeshRequest &sourceRequest,
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
snap_kind_priority(SoBRLSnapAction::SnapKind kind)
{
    switch (kind) {
	case SoBRLSnapAction::ENDPOINT:
	    return 0;
	case SoBRLSnapAction::MIDPOINT:
	    return 1;
	case SoBRLSnapAction::CENTER:
	    return 2;
	case SoBRLSnapAction::FACE_NEAREST:
	    return 3;
	case SoBRLSnapAction::LINE_NEAREST:
	    return 4;
	case SoBRLSnapAction::CONSTRUCTION_PLANE:
	    return 5;
	case SoBRLSnapAction::NONE:
	default:
	    return 6;
    }
}

static float
distance_between(const SbVec3f &a, const SbVec3f &b)
{
    return (a - b).length();
}

static float
snap_source_local_tolerance(const SbMatrix &worldToLocal, float tolerance)
{
    if (tolerance <= 0.0f)
	return tolerance;

    const SbVec3f axes[3] = {
	SbVec3f(1.0f, 0.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f),
	SbVec3f(0.0f, 0.0f, 1.0f)
    };
    float scale = 0.0f;
    for (int i = 0; i < 3; i++) {
	SbVec3f localAxis;
	worldToLocal.multDirMatrix(axes[i], localAxis);
	float len = localAxis.length();
	if (len > scale)
	    scale = len;
    }

    return scale > 0.0f ? tolerance * scale : tolerance;
}

static SbVec3f
snap_closest_point_on_segment(const SbVec3f &query, const SbVec3f &a, const SbVec3f &b)
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
snap_closest_point_on_triangle(const SbVec3f &query,
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

static SbBool
closest_point_on_plane(const SbVec3f &query,
	const SbVec3f &origin,
	const SbVec3f &normal,
	SbVec3f &candidate)
{
    SbVec3f n = normal;
    float len = n.length();
    if (len <= 0.0f)
	return FALSE;
    n /= len;
    candidate = query - n * ((query - origin).dot(n));
    return TRUE;
}

SoBRLSnapAction::SoBRLSnapAction(void) :
    queryPoint(0.0f, 0.0f, 0.0f),
    candidatePoint(0.0f, 0.0f, 0.0f),
    candidatePath(""),
    enabledKinds(ALL_KINDS),
    tolerance(0.25f),
    bestDistance(FLT_MAX),
    candidatePrimitiveIndex(-1),
    candidateKind(NONE),
    constructionPlaneOrigin(0.0f, 0.0f, 0.0f),
    constructionPlaneNormal(0.0f, 0.0f, 1.0f),
    constructionPlanePath("construction::plane"),
    selectionFilter(ALL_GEOMETRY),
    coordinateSpace(WORLD_SPACE),
    priorityPolicy(NEAREST_DISTANCE),
    geometryPolicy(SoBRLSnapAction::DISPLAY_LEVEL),
    skippedLodDisplayMeshCount(0),
    foundCandidate(FALSE),
    constructionPlaneEnabled(FALSE)
{
    SO_ACTION_CONSTRUCTOR(SoBRLSnapAction);
}

SoBRLSnapAction::~SoBRLSnapAction(void)
{
}

void
SoBRLSnapAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLSnapAction, SoAction);
    SO_ENABLE(SoBRLSnapAction, SoModelMatrixElement);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLSnapAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLVListShape, SoBRLSnapAction::vlistShapeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape, SoBRLSnapAction::meshShapeAction);
}

void
SoBRLSnapAction::setQueryPoint(const SbVec3f &point)
{
    this->queryPoint = point;
}

void
SoBRLSnapAction::setTolerance(float tol)
{
    this->tolerance = tol;
}

void
SoBRLSnapAction::setEnabledKinds(uint32_t kinds)
{
    this->enabledKinds = kinds;
}

uint32_t
SoBRLSnapAction::getEnabledKinds(void) const
{
    return this->enabledKinds;
}

void
SoBRLSnapAction::setSelectionFilter(SelectionFilter filter)
{
    this->selectionFilter = filter;
}

SoBRLSnapAction::SelectionFilter
SoBRLSnapAction::getSelectionFilter(void) const
{
    return this->selectionFilter;
}

void
SoBRLSnapAction::setCoordinateSpace(CoordinateSpace space)
{
    this->coordinateSpace = space;
}

SoBRLSnapAction::CoordinateSpace
SoBRLSnapAction::getCoordinateSpace(void) const
{
    return this->coordinateSpace;
}

void
SoBRLSnapAction::setPriorityPolicy(PriorityPolicy policy)
{
    this->priorityPolicy = policy;
}

SoBRLSnapAction::PriorityPolicy
SoBRLSnapAction::getPriorityPolicy(void) const
{
    return this->priorityPolicy;
}

void
SoBRLSnapAction::setGeometryPolicy(GeometryPolicy policy)
{
    this->geometryPolicy = (policy == SoBRLSnapAction::DISPLAY_LEVEL) ?
	SoBRLSnapAction::DISPLAY_LEVEL : SoBRLSnapAction::FULL_DETAIL;
}

SoBRLSnapAction::GeometryPolicy
SoBRLSnapAction::getGeometryPolicy(void) const
{
    return this->geometryPolicy;
}

unsigned int
SoBRLSnapAction::getSkippedLodDisplayMeshCount(void) const
{
    return this->skippedLodDisplayMeshCount;
}

int
SoBRLSnapAction::getSourceBackedFullDetailRequestCount(void) const
{
    return static_cast<int>(this->sourceBackedFullDetailRequests.size());
}

const BRLObolSourceMeshRequest &
SoBRLSnapAction::getSourceBackedFullDetailRequest(int index) const
{
    return this->sourceBackedFullDetailRequests.at(static_cast<size_t>(index));
}

SbBool
SoBRLSnapAction::makeSourceBackedFullDetailLodRequest(int index,
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
SoBRLSnapAction::consumeSourceBackedFullDetailResult(
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result)
{
    if (!snap_source_full_detail_result_valid(sourceRequest, result))
	return FALSE;
    if (!this->selectionAllows(sourceRequest.selected))
	return TRUE;

    const SbVec3f query = this->queryPoint;
    SbBox3f centerBox;
    centerBox.makeEmpty();
    size_t faceCount = result.mesh.coordIndex.size() / 3;
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
	centerBox.extendBy(pointA);
	centerBox.extendBy(pointB);
	centerBox.extendBy(pointC);

	this->consider(FACE_NEAREST, sourceRequest.path, static_cast<int>(i),
		query, snap_closest_point_on_triangle(query, pointA, pointB,
		    pointC));
    }

    if (!centerBox.isEmpty())
	this->consider(CENTER, sourceRequest.path, -1, query,
		centerBox.getCenter());

    return TRUE;
}

int
SoBRLSnapAction::submitSourceBackedFullDetailRequests(
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
SoBRLSnapAction::consumeSourceBackedFullDetailResults(
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

void
SoBRLSnapAction::setConstructionPlane(const SbVec3f &origin, const SbVec3f &normal)
{
    this->setConstructionPlane(origin, normal, SbString("construction::plane"));
}

void
SoBRLSnapAction::setConstructionPlane(const SbVec3f &origin,
	const SbVec3f &normal,
	const SbString &path)
{
    this->constructionPlaneOrigin = origin;
    this->constructionPlaneNormal = normal;
    this->constructionPlanePath = path;
    this->constructionPlaneEnabled = TRUE;
}

void
SoBRLSnapAction::clearConstructionPlane(void)
{
    this->constructionPlaneEnabled = FALSE;
}

SbBool
SoBRLSnapAction::hasConstructionPlane(void) const
{
    return this->constructionPlaneEnabled;
}

SbBool
SoBRLSnapAction::hasCandidate(void) const
{
    return this->foundCandidate;
}

SoBRLSnapAction::SnapKind
SoBRLSnapAction::getKind(void) const
{
    return this->candidateKind;
}

const SbVec3f &
SoBRLSnapAction::getPoint(void) const
{
    return this->candidatePoint;
}

const SbString &
SoBRLSnapAction::getPath(void) const
{
    return this->candidatePath;
}

int
SoBRLSnapAction::getPrimitiveIndex(void) const
{
    return this->candidatePrimitiveIndex;
}

float
SoBRLSnapAction::getDistance(void) const
{
    return this->bestDistance;
}

void
SoBRLSnapAction::beginTraversal(SoNode *node)
{
    this->candidatePoint = SbVec3f(0.0f, 0.0f, 0.0f);
    this->candidatePath = "";
    this->bestDistance = FLT_MAX;
    this->candidatePrimitiveIndex = -1;
    this->candidateKind = NONE;
    this->foundCandidate = FALSE;
    this->skippedLodDisplayMeshCount = 0;
    this->sourceBackedFullDetailRequests.clear();
    this->traverse(node);
    if (this->constructionPlaneEnabled) {
	SbVec3f planePoint;
	if (closest_point_on_plane(this->queryPoint, this->constructionPlaneOrigin,
		this->constructionPlaneNormal, planePoint)) {
	    this->consider(CONSTRUCTION_PLANE, this->constructionPlanePath, -1,
		    this->queryPoint, planePoint);
	}
    }
}

void
SoBRLSnapAction::nodeAction(SoAction *action, SoNode *node)
{
    node->doAction(action);
}

void
SoBRLSnapAction::vlistShapeAction(SoAction *action, SoNode *node)
{
    SoBRLSnapAction *snapAction = static_cast<SoBRLSnapAction *>(action);
    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    if (!shape->visible.getValue() || !shape->selectable.getValue())
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    const SbVec3f query = snapAction->queryPoint;
    SbBox3f centerBox;
    centerBox.makeEmpty();

    for (int i = 0; i < shape->getSegmentCount(); i++) {
	SbVec3f a;
	SbVec3f b;
	if (!shape->getSegment(i, a, b))
	    continue;
	if (!snapAction->selectionAllows(shape->isPrimitiveSelected(i)))
	    continue;

	SbVec3f pointA = snapAction->pointForCoordinateSpace(localToWorld, a);
	SbVec3f pointB = snapAction->pointForCoordinateSpace(localToWorld, b);
	centerBox.extendBy(pointA);
	centerBox.extendBy(pointB);

	snapAction->consider(ENDPOINT, shape->sourcePath.getValue(), i, query, pointA);
	snapAction->consider(ENDPOINT, shape->sourcePath.getValue(), i, query, pointB);
	snapAction->consider(LINE_NEAREST, shape->sourcePath.getValue(), i, query,
		snap_closest_point_on_segment(query, pointA, pointB));
	snapAction->consider(MIDPOINT, shape->sourcePath.getValue(), i, query,
		(pointA + pointB) * 0.5f);
    }

    for (int i = 0; i < shape->getPointPrimitiveCount(); i++) {
	SbVec3f point;
	int primitiveIndex = -1;
	if (!shape->getPointPrimitive(i, primitiveIndex, point))
	    continue;
	if (!snapAction->selectionAllows(shape->isPrimitiveSelected(primitiveIndex)))
	    continue;

	SbVec3f worldPoint = snapAction->pointForCoordinateSpace(localToWorld, point);
	centerBox.extendBy(worldPoint);
	snapAction->consider(ENDPOINT, shape->sourcePath.getValue(),
		primitiveIndex, query, worldPoint);
    }

    if (!centerBox.isEmpty())
	snapAction->consider(CENTER, shape->sourcePath.getValue(), -1, query,
		centerBox.getCenter());
}

void
SoBRLSnapAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLSnapAction *snapAction = static_cast<SoBRLSnapAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    if (!shape->visible.getValue() || !shape->selectable.getValue())
	return;
    const SbBool useFullDetail =
	(snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	 shape->hasFullDetailMesh()) ? TRUE : FALSE;
    const SbBool useSourceBackedFullDetail =
	(snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	 shape->needsSourceBackedFullDetail()) ? TRUE : FALSE;
    if (snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	    shape->isLodDisplayActive())
	snapAction->skippedLodDisplayMeshCount++;
    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    if (useSourceBackedFullDetail) {
	snapAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
	return;
    }
    const SbVec3f query = snapAction->queryPoint;
    SbBox3f centerBox;
    centerBox.makeEmpty();

    int triangleCount = useFullDetail ?
	shape->getFullDetailTriangleCount() : shape->getTriangleCount();
    for (int i = 0; i < triangleCount; i++) {
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
	if (!snapAction->selectionAllows(shape->isPrimitiveSelected(i)))
	    continue;

	SbVec3f pointA = snapAction->pointForCoordinateSpace(localToWorld, a);
	SbVec3f pointB = snapAction->pointForCoordinateSpace(localToWorld, b);
	SbVec3f pointC = snapAction->pointForCoordinateSpace(localToWorld, c);
	centerBox.extendBy(pointA);
	centerBox.extendBy(pointB);
	centerBox.extendBy(pointC);

	snapAction->consider(FACE_NEAREST, shape->sourcePath.getValue(), i, query,
		snap_closest_point_on_triangle(query, pointA, pointB, pointC));
    }

    if (!centerBox.isEmpty())
	snapAction->consider(CENTER, shape->sourcePath.getValue(), -1, query,
		centerBox.getCenter());
}

void
SoBRLSnapAction::appendSourceBackedFullDetailRequest(
	const SoBRLMeshShape *shape, const SbMatrix &localToWorld)
{
    if (!shape)
	return;

    BRLObolSourceMeshRequest request;
    if (shape->makeSourceMeshRequest(request)) {
	request.localToWorld = localToWorld;
	SbVec3f localQuery = this->queryPoint;
	float localTolerance = this->tolerance;
	if (this->coordinateSpace == WORLD_SPACE) {
	    SbMatrix worldToLocal = localToWorld.inverse();
	    worldToLocal.multVecMatrix(this->queryPoint, localQuery);
	    localTolerance = snap_source_local_tolerance(worldToLocal,
		    this->tolerance);
	}
	if (localTolerance >= 0.0f) {
	    SbVec3f delta(localTolerance, localTolerance, localTolerance);
	    request.queryBoundsValid = 1;
	    request.queryBounds.makeEmpty();
	    request.queryBounds.extendBy(localQuery - delta);
	    request.queryBounds.extendBy(localQuery + delta);
	    request.queryToleranceValid = 1;
	    request.queryTolerance = localTolerance;
	}
	this->sourceBackedFullDetailRequests.push_back(request);
    }
}

void
SoBRLSnapAction::consider(SnapKind kind, const SbString &path, int primitiveIndex,
	const SbVec3f &query, const SbVec3f &candidate)
{
    const float tieTolerance = 1.0e-6f;

    if (!(this->enabledKinds & static_cast<uint32_t>(kind)))
	return;

    float dist = distance_between(query, candidate);
    if (dist > this->tolerance)
	return;
    if (this->foundCandidate) {
	if (dist > this->bestDistance + tieTolerance)
	    return;
	if (fabsf(dist - this->bestDistance) <= tieTolerance) {
	    if (this->priorityPolicy != FEATURE_PRIORITY ||
		    snap_kind_priority(kind) >= snap_kind_priority(this->candidateKind))
		return;
	}
    }

    this->foundCandidate = TRUE;
    this->candidateKind = kind;
    this->candidatePath = path;
    this->candidatePrimitiveIndex = primitiveIndex;
    this->candidatePoint = candidate;
    this->bestDistance = dist;
}

SbVec3f
SoBRLSnapAction::pointForCoordinateSpace(const SbMatrix &localToWorld,
	const SbVec3f &localPoint) const
{
    if (this->coordinateSpace == PATH_LOCAL_SPACE)
	return localPoint;

    SbVec3f ret;
    localToWorld.multVecMatrix(localPoint, ret);
    return ret;
}

SbBool
SoBRLSnapAction::selectionAllows(SbBool selected) const
{
    if (this->selectionFilter == SELECTED_ONLY)
	return selected;
    if (this->selectionFilter == UNSELECTED_ONLY)
	return !selected;
    return TRUE;
}
