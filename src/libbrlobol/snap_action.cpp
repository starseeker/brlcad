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
#include "brlobol/view_lod.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SbBox.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoTransformation.h>

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
    if (sourceRequest.faceCount != 0) {
	const uint64_t resultFaceCount = static_cast<uint64_t>(faceCount);
	if (sourceRequest.queryBoundsValid && sourceRequest.queryToleranceValid) {
	    if (resultFaceCount > sourceRequest.faceCount)
		return FALSE;
	} else if (sourceRequest.faceCount != resultFaceCount) {
	    return FALSE;
	}
    }
    if (sourceRequest.pointCount != 0) {
	const uint64_t resultPointCount =
	    static_cast<uint64_t>(result.mesh.points.size());
	if (sourceRequest.queryBoundsValid && sourceRequest.queryToleranceValid) {
	    if (resultPointCount > sourceRequest.pointCount)
		return FALSE;
	    if (resultPointCount < sourceRequest.pointCount &&
		result.mesh.vertexIndex.size() != result.mesh.points.size())
		return FALSE;
	} else if (sourceRequest.pointCount != resultPointCount) {
	    return FALSE;
	}
    }

    return TRUE;
}

static int
snap_source_mesh_face_index(const BRLObolLodMeshPayload &mesh,
			    size_t faceSlot, size_t faceCount)
{
    if (mesh.faceIndex.size() == faceCount)
	return static_cast<int>(mesh.faceIndex[faceSlot]);
    return static_cast<int>(faceSlot);
}

static int
snap_source_mesh_vertex_index(const BRLObolLodMeshPayload &mesh,
			      int localIndex)
{
    if (localIndex < 0)
	return localIndex;
    if (mesh.vertexIndex.size() == mesh.points.size())
	return static_cast<int>(
		   mesh.vertexIndex[static_cast<size_t>(localIndex)]);
    return localIndex;
}

static int
snap_kind_priority(SoBRLSnapAction::SnapKind kind)
{
    switch (kind) {
	case SoBRLSnapAction::VERTEX:
	    return 0;
	case SoBRLSnapAction::ENDPOINT:
	    return 1;
	case SoBRLSnapAction::MIDPOINT:
	    return 2;
	case SoBRLSnapAction::EDGE_NEAREST:
	    return 3;
	case SoBRLSnapAction::FACE_NEAREST:
	    return 4;
	case SoBRLSnapAction::CENTER:
	    return 5;
	case SoBRLSnapAction::LINE_NEAREST:
	    return 6;
	case SoBRLSnapAction::CONSTRUCTION_PLANE:
	    return 7;
	case SoBRLSnapAction::NONE:
	default:
	    return 8;
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
    candidateEditIntentId(""),
    candidateEditIntentRole(""),
    enabledKinds(ALL_KINDS),
    tolerance(0.25f),
    bestDistance(FLT_MAX),
    candidatePrimitiveIndex(-1),
    candidateVertexIndex(-1),
    candidateEdgeSlot(-1),
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
    this->candidateEdgeVertexIndex[0] = -1;
    this->candidateEdgeVertexIndex[1] = -1;
}

SoBRLSnapAction::~SoBRLSnapAction(void)
{
}

void
SoBRLSnapAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLSnapAction, SoAction);
    SO_ENABLE(SoBRLSnapAction, SoModelMatrixElement);
    SO_ENABLE(SoBRLSnapAction, SoBRLViewLodElement);
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
	const int sourceFaceIndex =
	    snap_source_mesh_face_index(result.mesh, i, faceCount);
	const int sourceIa = snap_source_mesh_vertex_index(result.mesh, ia);
	const int sourceIb = snap_source_mesh_vertex_index(result.mesh, ib);
	const int sourceIc = snap_source_mesh_vertex_index(result.mesh, ic);
	centerBox.extendBy(pointA);
	centerBox.extendBy(pointB);
	centerBox.extendBy(pointC);

	this->consider(VERTEX, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query, pointA, sourceIa);
	this->consider(VERTEX, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query, pointB, sourceIb);
	this->consider(VERTEX, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query, pointC, sourceIc);
	this->consider(EDGE_NEAREST, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query,
		       snap_closest_point_on_segment(query, pointA, pointB),
		       -1, 0, sourceIa, sourceIb);
	this->consider(EDGE_NEAREST, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query,
		       snap_closest_point_on_segment(query, pointB, pointC),
		       -1, 1, sourceIb, sourceIc);
	this->consider(EDGE_NEAREST, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query,
		       snap_closest_point_on_segment(query, pointC, pointA),
		       -1, 2, sourceIc, sourceIa);
	this->consider(FACE_NEAREST, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       sourceFaceIndex, query,
		       snap_closest_point_on_triangle(query, pointA, pointB, pointC));
    }

    if (!centerBox.isEmpty())
	this->consider(CENTER, sourceRequest.path,
		       sourceRequest.editIntentId, sourceRequest.editIntentRole,
		       -1, query, centerBox.getCenter());

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

const SbString &
SoBRLSnapAction::getEditIntentId(void) const
{
    return this->candidateEditIntentId;
}

const SbString &
SoBRLSnapAction::getEditIntentRole(void) const
{
    return this->candidateEditIntentRole;
}

int
SoBRLSnapAction::getPrimitiveIndex(void) const
{
    return this->candidatePrimitiveIndex;
}

int
SoBRLSnapAction::getVertexIndex(void) const
{
    return this->candidateVertexIndex;
}

int
SoBRLSnapAction::getEdgeSlot(void) const
{
    return this->candidateEdgeSlot;
}

int
SoBRLSnapAction::getEdgeVertexIndexA(void) const
{
    return this->candidateEdgeVertexIndex[0];
}

int
SoBRLSnapAction::getEdgeVertexIndexB(void) const
{
    return this->candidateEdgeVertexIndex[1];
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
    this->candidateEditIntentId = "";
    this->candidateEditIntentRole = "";
    this->bestDistance = FLT_MAX;
    this->candidatePrimitiveIndex = -1;
    this->candidateVertexIndex = -1;
    this->candidateEdgeSlot = -1;
    this->candidateEdgeVertexIndex[0] = -1;
    this->candidateEdgeVertexIndex[1] = -1;
    this->candidateKind = NONE;
    this->foundCandidate = FALSE;
    this->skippedLodDisplayMeshCount = 0;
    this->sourceBackedFullDetailRequests.clear();
    this->traverse(node);
    if (this->constructionPlaneEnabled) {
	SbVec3f planePoint;
	if (closest_point_on_plane(this->queryPoint, this->constructionPlaneOrigin,
				   this->constructionPlaneNormal, planePoint)) {
	    this->consider(CONSTRUCTION_PLANE, this->constructionPlanePath,
			   SbString(""), SbString(""), -1, this->queryPoint, planePoint);
	}
    }
}

void
SoBRLSnapAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()) ||
	node->isOfType(SoTransformation::getClassTypeId()))
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
    const SbString &sourcePath = shape->sourcePath.getValue();
    const SbString &editIntentId = shape->editIntentId.getValue();
    const SbString &editIntentRole = shape->editIntentRole.getValue();
    SbBox3f centerBox;
    centerBox.makeEmpty();

    SbVec3f last;
    SbBool haveLast = FALSE;
    int segmentIndex = 0;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    int n = geom->point.getNum();
    if (geom->command.getNum() < n)
	n = geom->command.getNum();

    for (int i = 0; i < n; i++) {
	if (geom->command[i] == SoBRLVListShape::MOVE) {
	    last = geom->point[i];
	    haveLast = TRUE;
	    continue;
	}

	if (geom->command[i] == SoBRLVListShape::DRAW) {
	    if (!haveLast) {
		last = geom->point[i];
		haveLast = TRUE;
		continue;
	    }

	    const int currentSegment = segmentIndex++;
	    if (!snapAction->selectionAllows(shape->isPrimitiveSelected(currentSegment))) {
		last = geom->point[i];
		continue;
	    }

	    SbVec3f pointA = snapAction->pointForCoordinateSpace(localToWorld, last);
	    SbVec3f pointB = snapAction->pointForCoordinateSpace(localToWorld, geom->point[i]);
	    centerBox.extendBy(pointA);
	    centerBox.extendBy(pointB);

	    snapAction->consider(ENDPOINT, sourcePath, editIntentId,
				 editIntentRole, currentSegment, query, pointA);
	    snapAction->consider(ENDPOINT, sourcePath, editIntentId,
				 editIntentRole, currentSegment, query, pointB);
	    snapAction->consider(LINE_NEAREST, sourcePath, editIntentId,
				 editIntentRole, currentSegment, query,
				 snap_closest_point_on_segment(query, pointA, pointB));
	    snapAction->consider(MIDPOINT, sourcePath, editIntentId,
				 editIntentRole, currentSegment, query,
				 (pointA + pointB) * 0.5f);
	    last = geom->point[i];
	    continue;
	}

	if (geom->command[i] != SoBRLVListShape::POINT)
	    continue;

	if (!snapAction->selectionAllows(shape->isPrimitiveSelected(i)))
	    continue;

	SbVec3f worldPoint = snapAction->pointForCoordinateSpace(localToWorld, geom->point[i]);
	centerBox.extendBy(worldPoint);
	snapAction->consider(ENDPOINT, sourcePath, editIntentId,
			     editIntentRole, i, query, worldPoint);
    }

    if (!centerBox.isEmpty())
	snapAction->consider(CENTER, sourcePath, editIntentId, editIntentRole,
			     -1, query, centerBox.getCenter());
}

void
SoBRLSnapAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLSnapAction *snapAction = static_cast<SoBRLSnapAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    if (!shape->visible.getValue() || !shape->selectable.getValue())
	return;
    const BRLObolViewLodState::MeshPayload *viewPayload =
	brlobol_view_lod_mesh_for_action(action, shape);
    const SbBool useFullDetail =
	(snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	 shape->hasFullDetailMesh()) ? TRUE : FALSE;
    const SbBool useSourceBackedFullDetail =
	(snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	 shape->needsSourceBackedFullDetail()) ? TRUE : FALSE;
    if (snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	(shape->isLodDisplayActive() || viewPayload))
	snapAction->skippedLodDisplayMeshCount++;
    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    if (useSourceBackedFullDetail) {
	snapAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
	return;
    }
    const SbVec3f query = snapAction->queryPoint;
    const SbString &sourcePath = shape->sourcePath.getValue();
    const SbString &editIntentId = shape->editIntentId.getValue();
    const SbString &editIntentRole = shape->editIntentRole.getValue();
    SbBox3f centerBox;
    centerBox.makeEmpty();

    const SbBool useViewPayload =
	(!useFullDetail &&
	 snapAction->geometryPolicy == SoBRLSnapAction::DISPLAY_LEVEL &&
	 viewPayload) ? TRUE : FALSE;
    int triangleCount = useFullDetail ?
			shape->getFullDetailTriangleCount() :
			(useViewPayload ? viewPayload->getTriangleCount() :
			 shape->getTriangleCount());
    for (int i = 0; i < triangleCount; i++) {
	int ia = -1;
	int ib = -1;
	int ic = -1;
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (useFullDetail) {
	    if (!shape->getFullDetailTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    if (!shape->getFullDetailTriangle(i, a, b, c))
		continue;
	} else if (useViewPayload) {
	    if (!viewPayload->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    if (!viewPayload->getTriangle(i, a, b, c))
		continue;
	} else {
	    if (!shape->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
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

	snapAction->consider(VERTEX, sourcePath, editIntentId, editIntentRole,
			     i, query, pointA, ia);
	snapAction->consider(VERTEX, sourcePath, editIntentId, editIntentRole,
			     i, query, pointB, ib);
	snapAction->consider(VERTEX, sourcePath, editIntentId, editIntentRole,
			     i, query, pointC, ic);
	snapAction->consider(EDGE_NEAREST, sourcePath, editIntentId,
			     editIntentRole, i, query,
			     snap_closest_point_on_segment(query, pointA, pointB),
			     -1, 0, ia, ib);
	snapAction->consider(EDGE_NEAREST, sourcePath, editIntentId,
			     editIntentRole, i, query,
			     snap_closest_point_on_segment(query, pointB, pointC),
			     -1, 1, ib, ic);
	snapAction->consider(EDGE_NEAREST, sourcePath, editIntentId,
			     editIntentRole, i, query,
			     snap_closest_point_on_segment(query, pointC, pointA),
			     -1, 2, ic, ia);
	snapAction->consider(FACE_NEAREST, sourcePath, editIntentId,
			     editIntentRole, i, query,
			     snap_closest_point_on_triangle(query, pointA, pointB, pointC));
    }

    if (!centerBox.isEmpty())
	snapAction->consider(CENTER, sourcePath, editIntentId, editIntentRole,
			     -1, query, centerBox.getCenter());
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
SoBRLSnapAction::consider(SnapKind kind, const SbString &path,
			  const SbString &editIntentId,
			  const SbString &editIntentRole,
			  int primitiveIndex,
			  const SbVec3f &query, const SbVec3f &candidate,
			  int vertexIndex,
			  int edgeSlot,
			  int edgeVertexIndexA,
			  int edgeVertexIndexB)
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
    this->candidateEditIntentId = editIntentId;
    this->candidateEditIntentRole = editIntentRole;
    this->candidatePrimitiveIndex = primitiveIndex;
    this->candidateVertexIndex = vertexIndex;
    this->candidateEdgeSlot = edgeSlot;
    this->candidateEdgeVertexIndex[0] = edgeVertexIndexA;
    this->candidateEdgeVertexIndex[1] = edgeVertexIndexB;
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
