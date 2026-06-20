/*                    S N A P _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/mesh_shape.h"
#include "brlobol/snap_action.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SbBox.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

#include <float.h>
#include <math.h>

SO_ACTION_SOURCE(SoBRLSnapAction);

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
	 shape->isLodDisplayActive() && shape->hasFullDetailMesh()) ?
	TRUE : FALSE;
    if (snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	    shape->isLodDisplayActive())
	snapAction->skippedLodDisplayMeshCount++;
    if (snapAction->geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	    shape->isLodDisplayActive() && !useFullDetail)
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
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
