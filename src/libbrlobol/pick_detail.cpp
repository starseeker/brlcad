/*                   P I C K _ D E T A I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_realization.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/pick_detail.h"

#include "raytrace.h"

#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoNode.h>

#include <float.h>
#include <math.h>
#include <vector>

SO_DETAIL_SOURCE(SoBRLPickDetail);
SO_ACTION_SOURCE(SoBRLSourceMeshPickAction);

BRLObolSourceMeshPickResult::BRLObolSourceMeshPickResult(void)
{
    clear();
}

void
BRLObolSourceMeshPickResult::clear(void)
{
    detail.clear();
    point = SbVec3f(0.0f, 0.0f, 0.0f);
    distance = FLT_MAX;
    hit = FALSE;
}

BRLObolRtPickResult::BRLObolRtPickResult(void)
{
    clear();
}

void
BRLObolRtPickResult::clear(void)
{
    detail.clear();
    point = SbVec3f(0.0f, 0.0f, 0.0f);
    normal = SbVec3f(0.0f, 0.0f, 0.0f);
    distance = FLT_MAX;
    hit = FALSE;
}

static float
pick_source_distance_squared_to_segment(const SbVec3f &p, const SbVec3f &a,
	const SbVec3f &b)
{
    SbVec3f ab = b - a;
    float denom = ab.sqrLength();
    if (denom <= 0.0f)
	return (p - a).sqrLength();

    float t = (p - a).dot(ab) / denom;
    if (t < 0.0f)
	t = 0.0f;
    else if (t > 1.0f)
	t = 1.0f;

    return (p - (a + ab * t)).sqrLength();
}

static int
pick_source_nearest_face_vertex_slot(const SbVec3f &hit,
	const SbVec3f vertices[3])
{
    int nearest = 0;
    float nearestDist = (hit - vertices[0]).sqrLength();
    for (int i = 1; i < 3; i++) {
	float dist = (hit - vertices[i]).sqrLength();
	if (dist < nearestDist) {
	    nearest = i;
	    nearestDist = dist;
	}
    }
    return nearest;
}

static int
pick_source_nearest_face_edge_slot(const SbVec3f &hit,
	const SbVec3f vertices[3])
{
    static const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    int nearest = 0;
    float nearestDist = pick_source_distance_squared_to_segment(hit,
	    vertices[edges[0][0]], vertices[edges[0][1]]);
    for (int i = 1; i < 3; i++) {
	float dist = pick_source_distance_squared_to_segment(hit,
		vertices[edges[i][0]], vertices[edges[i][1]]);
	if (dist < nearestDist) {
	    nearest = i;
	    nearestDist = dist;
	}
    }
    return nearest;
}

static SbBool
pick_source_ray_triangle(float &t, float &u, float &v,
	const SbVec3f &origin, const SbVec3f &direction,
	const SbVec3f &a, const SbVec3f &b, const SbVec3f &c)
{
    const float epsilon = 1.0e-7f;
    SbVec3f ab = b - a;
    SbVec3f ac = c - a;
    SbVec3f pvec = direction.cross(ac);
    float det = ab.dot(pvec);
    if (fabsf(det) < epsilon)
	return FALSE;

    float invDet = 1.0f / det;
    SbVec3f tvec = origin - a;
    u = tvec.dot(pvec) * invDet;
    if (u < 0.0f || u > 1.0f)
	return FALSE;

    SbVec3f qvec = tvec.cross(ab);
    v = direction.dot(qvec) * invDet;
    if (v < 0.0f || u + v > 1.0f)
	return FALSE;

    t = ac.dot(qvec) * invDet;
    return t >= 0.0f ? TRUE : FALSE;
}

static SbBool
pick_source_full_detail_result_valid(
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result)
{
    if (result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	    result.resultKind != BRLOBOL_LOD_RESULT_FULL_DETAIL ||
	    !result.mesh.isValid())
	return FALSE;

    size_t faceCount = result.mesh.coordIndex.size() / 3;
    if (sourceRequest.faceCount != 0) {
	const uint64_t resultFaceCount = static_cast<uint64_t>(faceCount);
	if (sourceRequest.queryRayValid) {
	    if (resultFaceCount > sourceRequest.faceCount)
		return FALSE;
	} else if (sourceRequest.faceCount != resultFaceCount) {
	    return FALSE;
	}
    }
    if (sourceRequest.pointCount != 0) {
	const uint64_t resultPointCount =
	    static_cast<uint64_t>(result.mesh.points.size());
	if (sourceRequest.queryRayValid) {
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
pick_source_mesh_face_index(const BRLObolLodMeshPayload &mesh,
	size_t faceSlot, size_t faceCount)
{
    if (mesh.faceIndex.size() == faceCount)
	return static_cast<int>(mesh.faceIndex[faceSlot]);
    return static_cast<int>(faceSlot);
}

static int
pick_source_mesh_vertex_index(const BRLObolLodMeshPayload &mesh,
	int localIndex)
{
    if (localIndex < 0)
	return localIndex;
    if (mesh.vertexIndex.size() == mesh.points.size())
	return static_cast<int>(
		mesh.vertexIndex[static_cast<size_t>(localIndex)]);
    return localIndex;
}

static SbString
pick_path_string(const char *name)
{
    if (!name || !name[0])
	return "";
    if (name[0] == '/')
	return SbString(name);

    SbString path("/");
    path += name;
    return path;
}

struct brlobol_rt_pick_state {
    BRLObolRtPickResult pick;
};

static int
brlobol_pick_rt_hit(struct application *ap,
	struct partition *PartHeadp,
	struct seg *UNUSED(segp))
{
    if (!ap || !ap->a_uptr || !PartHeadp)
	return 0;

    brlobol_rt_pick_state *state =
	static_cast<brlobol_rt_pick_state *>(ap->a_uptr);

    for (struct partition *pp = PartHeadp->pt_forw;
	    pp != PartHeadp; pp = pp->pt_forw) {
	if (!pp || !pp->pt_inhit || pp->pt_inhit->hit_dist < 0.0 ||
		pp->pt_inhit->hit_dist >= state->pick.distance)
	    continue;

	struct soltab *stp = pp->pt_inseg ? pp->pt_inseg->seg_stp : NULL;
	point_t hitPoint = VINIT_ZERO;
	vect_t hitNormal = VINIT_ZERO;
	VJOIN1(hitPoint, ap->a_ray.r_pt, pp->pt_inhit->hit_dist,
		ap->a_ray.r_dir);
	if (stp)
	    RT_HIT_NORMAL(hitNormal, pp->pt_inhit, stp, &(ap->a_ray),
		    pp->pt_inflip);

	state->pick.clear();
	state->pick.hit = TRUE;
	state->pick.distance = static_cast<float>(pp->pt_inhit->hit_dist);
	state->pick.point = SbVec3f(static_cast<float>(hitPoint[X]),
		static_cast<float>(hitPoint[Y]),
		static_cast<float>(hitPoint[Z]));
	state->pick.normal = SbVec3f(static_cast<float>(hitNormal[X]),
		static_cast<float>(hitNormal[Y]),
		static_cast<float>(hitNormal[Z]));

	const char *regionName = pp->pt_regionp ?
	    pp->pt_regionp->reg_name : NULL;
	const char *primitiveName = stp && stp->st_dp ?
	    stp->st_dp->d_namep : NULL;
	const char *typeName = stp && stp->st_meth ?
	    stp->st_meth->ft_label : NULL;
	state->pick.detail.setPath(pick_path_string(
		regionName ? regionName : primitiveName));
	state->pick.detail.setSourceName(primitiveName ? primitiveName : "");
	state->pick.detail.setSourceType(typeName ? typeName : "");
	state->pick.detail.setPrimitive(SoBRLPickDetail::IMPLICIT_SOLID,
		pp->pt_inhit->hit_surfno);
	state->pick.detail.setModelPoint(state->pick.point);
	if (pp->pt_regionp) {
	    state->pick.detail.setRegionId(pp->pt_regionp->reg_regionid);
	    state->pick.detail.setAirCode(pp->pt_regionp->reg_aircode);
	    state->pick.detail.setMaterialId(pp->pt_regionp->reg_gmater);
	    state->pick.detail.setLos(pp->pt_regionp->reg_los);
	    if (pp->pt_regionp->reg_mater.ma_color_valid) {
		state->pick.detail.setMaterialColor(TRUE,
			SbColor(
			    static_cast<float>(pp->pt_regionp->reg_mater.ma_color[0]),
			    static_cast<float>(pp->pt_regionp->reg_mater.ma_color[1]),
			    static_cast<float>(pp->pt_regionp->reg_mater.ma_color[2])));
	    }
	    if (pp->pt_regionp->reg_mater.ma_shader)
		state->pick.detail.setMaterialShader(
			pp->pt_regionp->reg_mater.ma_shader);
	}
    }

    return state->pick.hit ? 1 : 0;
}

static int
brlobol_pick_rt_miss(struct application *UNUSED(ap))
{
    return 0;
}

static SbString
pick_rt_normalized_object_path(const SbString &objectPath)
{
    const char *name = objectPath.getString();
    while (name && name[0] == '/')
	name++;
    return name && name[0] ? SbString(name) : SbString("");
}

static SbBool
pick_rt_same_object_paths(const std::vector<SbString> &a,
	const std::vector<SbString> &b)
{
    if (a.size() != b.size())
	return FALSE;
    for (size_t i = 0; i < a.size(); i++) {
	if (strcmp(a[i].getString(), b[i].getString()) != 0)
	    return FALSE;
    }
    return TRUE;
}

BRLObolRtPickCache::BRLObolRtPickCache(void) :
    rtip(NULL),
    database(NULL),
    ready(FALSE)
{
}

BRLObolRtPickCache::~BRLObolRtPickCache(void)
{
    this->clear();
}

void
BRLObolRtPickCache::clear(void)
{
    if (this->rtip)
	rt_i_destroy(this->rtip);
    this->rtip = NULL;
    this->database = NULL;
    this->objectPaths.clear();
    this->ready = FALSE;
}

SbBool
BRLObolRtPickCache::prepare(struct db_i *dbip,
	const std::vector<SbString> &paths)
{
    if (!dbip || paths.empty()) {
	this->clear();
	return FALSE;
    }

    std::vector<SbString> normalizedPaths;
    normalizedPaths.reserve(paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
	SbString normalized = pick_rt_normalized_object_path(paths[i]);
	const char *name = normalized.getString();
	if (name && name[0] &&
		db_lookup(dbip, name, LOOKUP_QUIET) != RT_DIR_NULL)
	    normalizedPaths.push_back(normalized);
    }
    if (normalizedPaths.empty()) {
	this->clear();
	return FALSE;
    }

    if (this->ready && this->rtip && this->database == dbip &&
	    pick_rt_same_object_paths(this->objectPaths, normalizedPaths))
	return TRUE;

    std::vector<const char *> names;
    names.reserve(normalizedPaths.size());
    for (size_t i = 0; i < normalizedPaths.size(); i++)
	names.push_back(normalizedPaths[i].getString());

    struct rt_i *newRtip = rt_i_create(dbip);
    if (!newRtip) {
	this->clear();
	return FALSE;
    }
    newRtip->useair = 1;
    if (rt_gettrees(newRtip, static_cast<int>(names.size()), &names[0],
	    1) < 0) {
	rt_i_destroy(newRtip);
	this->clear();
	return FALSE;
    }
    newRtip->rti_hasty_prep = 1;
    rt_prep(newRtip);

    this->clear();
    this->rtip = newRtip;
    this->database = dbip;
    this->objectPaths = normalizedPaths;
    this->ready = TRUE;
    return TRUE;
}

SbBool
BRLObolRtPickCache::isReady(void) const
{
    return this->ready && this->rtip ? TRUE : FALSE;
}

int
BRLObolRtPickCache::getObjectPathCount(void) const
{
    return static_cast<int>(this->objectPaths.size());
}

SbBool
BRLObolRtPickCache::pickRay(BRLObolRtPickResult &pick,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection) const
{
    pick.clear();
    if (!this->isReady())
	return FALSE;

    SbVec3f direction = rayDirection;
    if (direction.length() <= 0.0f)
	return FALSE;
    direction.normalize();

    brlobol_rt_pick_state state;
    struct application ap;
    RT_APPLICATION_INIT(&ap);
    ap.a_magic = RT_AP_MAGIC;
    ap.a_ray.magic = RT_RAY_MAGIC;
    ap.a_hit = brlobol_pick_rt_hit;
    ap.a_miss = brlobol_pick_rt_miss;
    ap.a_overlap = NULL;
    ap.a_onehit = 0;
    ap.a_rt_i = this->rtip;
    ap.a_uptr = &state;
    ap.a_purpose = "BRLObolRtPickCache::pickRay";
    VSET(ap.a_ray.r_pt, rayOrigin[0], rayOrigin[1], rayOrigin[2]);
    VSET(ap.a_ray.r_dir, direction[0], direction[1], direction[2]);
    (void)rt_shootray(&ap);

    pick = state.pick;
    return pick.hit;
}

SbBool
brlobol_pick_rt_ray(BRLObolRtPickResult &pick,
	struct db_i *dbip,
	const std::vector<SbString> &objectPaths,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection)
{
    BRLObolRtPickCache cache;
    if (!cache.prepare(dbip, objectPaths)) {
	pick.clear();
	return FALSE;
    }

    return cache.pickRay(pick, rayOrigin, rayDirection);
}

static void
pick_source_fill_detail(SoBRLPickDetail &detail,
	const BRLObolSourceMeshRequest &sourceRequest,
	int primitiveIndex,
	const int vertexIndices[3],
	const SbVec3f localVertices[3],
	const SbVec3f &localHit)
{
    static const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

    detail.setPath(sourceRequest.path);
    detail.setSourceName(sourceRequest.sourceName);
    detail.setSourceType(sourceRequest.sourceType);
    detail.setSourceId(sourceRequest.sourceId);
    detail.setRegionId(sourceRequest.regionId);
    detail.setAirCode(sourceRequest.airCode);
    detail.setMaterialId(sourceRequest.materialId);
    detail.setLos(sourceRequest.los);
    detail.setMaterialColor(sourceRequest.materialColorValid ? TRUE : FALSE,
	    sourceRequest.materialColor);
    detail.setMaterialShader(sourceRequest.materialShader);
    detail.setPrimitive(SoBRLPickDetail::FACE, primitiveIndex);
    detail.setEditIntent(sourceRequest.editIntentId,
	    sourceRequest.editIntentRole);
    detail.setFaceVertexIndices(vertexIndices[0], vertexIndices[1],
	    vertexIndices[2]);
    int vertexSlot = pick_source_nearest_face_vertex_slot(localHit,
	    localVertices);
    int edgeSlot = pick_source_nearest_face_edge_slot(localHit, localVertices);
    detail.setNearestFaceVertex(vertexSlot, vertexIndices[vertexSlot]);
    detail.setNearestFaceEdge(edgeSlot, vertexIndices[edges[edgeSlot][0]],
	    vertexIndices[edges[edgeSlot][1]]);
    detail.setModelPoint(localHit);
}

SbBool
brlobol_pick_source_full_detail_result(
	BRLObolSourceMeshPickResult &pick,
	const BRLObolSourceMeshRequest &sourceRequest,
	const BRLObolLodResult &result,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection)
{
    pick.clear();
    if (!pick_source_full_detail_result_valid(sourceRequest, result))
	return FALSE;

    SbVec3f direction = rayDirection;
    if (direction.length() <= 0.0f)
	return FALSE;
    direction.normalize();

    size_t faceCount = result.mesh.coordIndex.size() / 3;
    for (size_t i = 0; i < faceCount; i++) {
	int localVertexIndices[3] = {
	    result.mesh.coordIndex[i * 3],
	    result.mesh.coordIndex[i * 3 + 1],
	    result.mesh.coordIndex[i * 3 + 2]
	};
	if (localVertexIndices[0] < 0 || localVertexIndices[1] < 0 ||
		localVertexIndices[2] < 0 ||
		static_cast<size_t>(localVertexIndices[0]) >=
		result.mesh.points.size() ||
		static_cast<size_t>(localVertexIndices[1]) >=
		result.mesh.points.size() ||
		static_cast<size_t>(localVertexIndices[2]) >=
		result.mesh.points.size())
	    return FALSE;

	SbVec3f localVertices[3] = {
	    result.mesh.points[static_cast<size_t>(localVertexIndices[0])],
	    result.mesh.points[static_cast<size_t>(localVertexIndices[1])],
	    result.mesh.points[static_cast<size_t>(localVertexIndices[2])]
	};
	int sourceVertexIndices[3] = {
	    pick_source_mesh_vertex_index(result.mesh, localVertexIndices[0]),
	    pick_source_mesh_vertex_index(result.mesh, localVertexIndices[1]),
	    pick_source_mesh_vertex_index(result.mesh, localVertexIndices[2])
	};
	SbVec3f worldVertices[3];
	for (int j = 0; j < 3; j++)
	    sourceRequest.localToWorld.multVecMatrix(localVertices[j],
		    worldVertices[j]);

	float t = 0.0f;
	float u = 0.0f;
	float v = 0.0f;
	if (!pick_source_ray_triangle(t, u, v, rayOrigin, direction,
		worldVertices[0], worldVertices[1], worldVertices[2]) ||
		t >= pick.distance)
	    continue;

	SbVec3f localHit = localVertices[0] +
	    (localVertices[1] - localVertices[0]) * u +
	    (localVertices[2] - localVertices[0]) * v;
	SbVec3f worldHit = worldVertices[0] +
	    (worldVertices[1] - worldVertices[0]) * u +
	    (worldVertices[2] - worldVertices[0]) * v;

	pick.distance = t;
	pick.point = worldHit;
	pick.hit = TRUE;
	pick_source_fill_detail(pick.detail, sourceRequest,
		pick_source_mesh_face_index(result.mesh, i, faceCount),
		sourceVertexIndices, localVertices, localHit);
    }

    return pick.hit;
}

static SbBool
pick_action_ray_intersects_box(const SbVec3f &origin, const SbVec3f &direction,
	const SbBox3f &box)
{
    if (box.isEmpty() || direction.length() <= 0.0f)
	return FALSE;

    float tmin = 0.0f;
    float tmax = FLT_MAX;
    const SbVec3f minp = box.getMin();
    const SbVec3f maxp = box.getMax();

    for (int axis = 0; axis < 3; axis++) {
	float o = origin[axis];
	float d = direction[axis];
	float mn = minp[axis];
	float mx = maxp[axis];

	if (fabsf(d) < 1.0e-8f) {
	    if (o < mn || o > mx)
		return FALSE;
	    continue;
	}

	float inv = 1.0f / d;
	float t1 = (mn - o) * inv;
	float t2 = (mx - o) * inv;
	if (t1 > t2) {
	    float tmp = t1;
	    t1 = t2;
	    t2 = tmp;
	}
	if (t1 > tmin)
	    tmin = t1;
	if (t2 < tmax)
	    tmax = t2;
	if (tmax < tmin)
	    return FALSE;
    }

    return tmax >= 0.0f ? TRUE : FALSE;
}

SoBRLSourceMeshPickAction::SoBRLSourceMeshPickAction(void) :
    rayOrigin(0.0f, 0.0f, 0.0f),
    rayDirection(0.0f, 0.0f, -1.0f),
    visitedMeshCount(0)
{
    SO_ACTION_CONSTRUCTOR(SoBRLSourceMeshPickAction);
}

SoBRLSourceMeshPickAction::~SoBRLSourceMeshPickAction(void)
{
}

void
SoBRLSourceMeshPickAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLSourceMeshPickAction, SoAction);
    SO_ENABLE(SoBRLSourceMeshPickAction, SoModelMatrixElement);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLSourceMeshPickAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape,
	    SoBRLSourceMeshPickAction::meshShapeAction);
}

void
SoBRLSourceMeshPickAction::setRay(const SbVec3f &origin,
	const SbVec3f &direction)
{
    this->rayOrigin = origin;
    this->rayDirection = direction;
}

const SbVec3f &
SoBRLSourceMeshPickAction::getRayOrigin(void) const
{
    return this->rayOrigin;
}

const SbVec3f &
SoBRLSourceMeshPickAction::getRayDirection(void) const
{
    return this->rayDirection;
}

int
SoBRLSourceMeshPickAction::getVisitedMeshCount(void) const
{
    return this->visitedMeshCount;
}

int
SoBRLSourceMeshPickAction::getSourceBackedFullDetailRequestCount(void) const
{
    return static_cast<int>(this->sourceBackedFullDetailRequests.size());
}

const BRLObolSourceMeshRequest &
SoBRLSourceMeshPickAction::getSourceBackedFullDetailRequest(int index) const
{
    return this->sourceBackedFullDetailRequests.at(static_cast<size_t>(index));
}

SbBool
SoBRLSourceMeshPickAction::makeSourceBackedFullDetailLodRequest(int index,
	BRLObolLodRequest &request,
	const BRLObolLodRequest *templateRequest) const
{
    if (index < 0 ||
	    static_cast<size_t>(index) >=
	    this->sourceBackedFullDetailRequests.size())
	return FALSE;

    return brlobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    request,
	    this->sourceBackedFullDetailRequests[static_cast<size_t>(index)],
	    templateRequest);
}

int
SoBRLSourceMeshPickAction::submitSourceBackedFullDetailRequests(
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
SoBRLSourceMeshPickAction::consumeSourceBackedFullDetailResults(
	BRLObolSourceMeshPickResult &pick,
	const std::vector<BRLObolLodResult> &results,
	const BRLObolLodRequest *templateRequest) const
{
    pick.clear();
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

	    BRLObolSourceMeshPickResult candidate;
	    if (brlobol_pick_source_full_detail_result(candidate,
		    this->sourceBackedFullDetailRequests[i], results[j],
		    this->rayOrigin, this->rayDirection)) {
		used[j] = TRUE;
		consumed++;
		if (!pick.hit || candidate.distance < pick.distance)
		    pick = candidate;
	    }
	    break;
	}
    }

    return consumed;
}

void
SoBRLSourceMeshPickAction::beginTraversal(SoNode *node)
{
    this->resetResults();
    this->traverse(node);
}

void
SoBRLSourceMeshPickAction::nodeAction(SoAction *action, SoNode *node)
{
    node->doAction(action);
}

void
SoBRLSourceMeshPickAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLSourceMeshPickAction *pickAction =
	static_cast<SoBRLSourceMeshPickAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);

    pickAction->visitedMeshCount++;
    if (!shape->visible.getValue() || !shape->selectable.getValue())
	return;
    if (shape->getPickGeometryPolicy() != SoBRLMeshShape::PICK_FULL_DETAIL ||
	    !shape->needsSourceBackedFullDetail())
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    pickAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
}

void
SoBRLSourceMeshPickAction::resetResults(void)
{
    this->sourceBackedFullDetailRequests.clear();
    this->visitedMeshCount = 0;
}

void
SoBRLSourceMeshPickAction::appendSourceBackedFullDetailRequest(
	const SoBRLMeshShape *shape, const SbMatrix &localToWorld)
{
    if (!shape || this->rayDirection.length() <= 0.0f)
	return;

    BRLObolSourceMeshRequest request;
    if (!shape->makeSourceMeshRequest(request))
	return;

    request.localToWorld = localToWorld;
    SbMatrix worldToLocal = localToWorld.inverse();
    worldToLocal.multVecMatrix(this->rayOrigin, request.queryRayOrigin);
    worldToLocal.multDirMatrix(this->rayDirection, request.queryRayDirection);
    if (request.queryRayDirection.length() > 0.0f) {
	request.queryRayDirection.normalize();
	request.queryRayValid = 1;
    }
    if (this->rayIntersectsBounds(request))
	this->sourceBackedFullDetailRequests.push_back(request);
}

SbBool
SoBRLSourceMeshPickAction::rayIntersectsBounds(
	const BRLObolSourceMeshRequest &request) const
{
    if (request.bounds.isEmpty())
	return TRUE;

    SbBox3f worldBox;
    worldBox.makeEmpty();
    const SbVec3f minp = request.bounds.getMin();
    const SbVec3f maxp = request.bounds.getMax();
    for (int xi = 0; xi < 2; xi++) {
	for (int yi = 0; yi < 2; yi++) {
	    for (int zi = 0; zi < 2; zi++) {
		SbVec3f localPoint(xi ? maxp[0] : minp[0],
			yi ? maxp[1] : minp[1],
			zi ? maxp[2] : minp[2]);
		SbVec3f worldPoint;
		request.localToWorld.multVecMatrix(localPoint, worldPoint);
		worldBox.extendBy(worldPoint);
	    }
	}
    }

    return pick_action_ray_intersects_box(this->rayOrigin,
	    this->rayDirection, worldBox);
}

SoBRLPickDetail::SoBRLPickDetail(void)
{
    clear();
}

SoBRLPickDetail::SoBRLPickDetail(const SoBRLPickDetail &other) :
    SoDetail(other),
    dbpath(other.dbpath),
    sourceName(other.sourceName),
    sourceType(other.sourceType),
    materialShader(other.materialShader),
    editIntentId(other.editIntentId),
    editIntentRole(other.editIntentRole),
    modelPoint(other.modelPoint),
    materialColor(other.materialColor),
    sourceId(other.sourceId),
    regionId(other.regionId),
    airCode(other.airCode),
    materialId(other.materialId),
    los(other.los),
    materialColorValid(other.materialColorValid),
    primitiveKind(other.primitiveKind),
    primitiveIndex(other.primitiveIndex)
{
    this->faceVertexIndex[0] = other.faceVertexIndex[0];
    this->faceVertexIndex[1] = other.faceVertexIndex[1];
    this->faceVertexIndex[2] = other.faceVertexIndex[2];
    this->nearestFaceEdgeSlot = other.nearestFaceEdgeSlot;
    this->nearestFaceEdgeVertexIndex[0] = other.nearestFaceEdgeVertexIndex[0];
    this->nearestFaceEdgeVertexIndex[1] = other.nearestFaceEdgeVertexIndex[1];
    this->nearestFaceVertexSlot = other.nearestFaceVertexSlot;
    this->nearestFaceVertexIndex = other.nearestFaceVertexIndex;
}

SoBRLPickDetail &
SoBRLPickDetail::operator=(const SoBRLPickDetail &other)
{
    if (this == &other)
	return *this;

    this->dbpath = other.dbpath;
    this->sourceName = other.sourceName;
    this->sourceType = other.sourceType;
    this->materialShader = other.materialShader;
    this->editIntentId = other.editIntentId;
    this->editIntentRole = other.editIntentRole;
    this->modelPoint = other.modelPoint;
    this->materialColor = other.materialColor;
    this->sourceId = other.sourceId;
    this->regionId = other.regionId;
    this->airCode = other.airCode;
    this->materialId = other.materialId;
    this->los = other.los;
    this->materialColorValid = other.materialColorValid;
    this->primitiveKind = other.primitiveKind;
    this->primitiveIndex = other.primitiveIndex;
    this->faceVertexIndex[0] = other.faceVertexIndex[0];
    this->faceVertexIndex[1] = other.faceVertexIndex[1];
    this->faceVertexIndex[2] = other.faceVertexIndex[2];
    this->nearestFaceEdgeSlot = other.nearestFaceEdgeSlot;
    this->nearestFaceEdgeVertexIndex[0] =
	other.nearestFaceEdgeVertexIndex[0];
    this->nearestFaceEdgeVertexIndex[1] =
	other.nearestFaceEdgeVertexIndex[1];
    this->nearestFaceVertexSlot = other.nearestFaceVertexSlot;
    this->nearestFaceVertexIndex = other.nearestFaceVertexIndex;

    return *this;
}

SoBRLPickDetail::~SoBRLPickDetail(void)
{
}

void
SoBRLPickDetail::initClass(void)
{
    SO_DETAIL_INIT_CLASS(SoBRLPickDetail, SoDetail);
}

SoDetail *
SoBRLPickDetail::copy(void) const
{
    return new SoBRLPickDetail(*this);
}

void
SoBRLPickDetail::clear(void)
{
    this->dbpath = "";
    this->sourceName = "";
    this->sourceType = "";
    this->materialShader = "";
    this->editIntentId = "";
    this->editIntentRole = "";
    this->modelPoint = SbVec3f(0.0f, 0.0f, 0.0f);
    this->materialColor = SbColor(1.0f, 1.0f, 1.0f);
    this->sourceId = 0;
    this->regionId = 0;
    this->airCode = 0;
    this->materialId = 0;
    this->los = 0;
    this->materialColorValid = FALSE;
    this->primitiveKind = UNKNOWN;
    this->primitiveIndex = -1;
    this->faceVertexIndex[0] = -1;
    this->faceVertexIndex[1] = -1;
    this->faceVertexIndex[2] = -1;
    this->nearestFaceEdgeSlot = -1;
    this->nearestFaceEdgeVertexIndex[0] = -1;
    this->nearestFaceEdgeVertexIndex[1] = -1;
    this->nearestFaceVertexSlot = -1;
    this->nearestFaceVertexIndex = -1;
}

void
SoBRLPickDetail::setPath(const SbString &path)
{
    this->dbpath = path;
}

const SbString &
SoBRLPickDetail::getPath(void) const
{
    return this->dbpath;
}

void
SoBRLPickDetail::setSourceName(const SbString &name)
{
    this->sourceName = name;
}

const SbString &
SoBRLPickDetail::getSourceName(void) const
{
    return this->sourceName;
}

void
SoBRLPickDetail::setSourceType(const SbString &type)
{
    this->sourceType = type;
}

const SbString &
SoBRLPickDetail::getSourceType(void) const
{
    return this->sourceType;
}

void
SoBRLPickDetail::setSourceId(uint32_t id)
{
    this->sourceId = id;
}

uint32_t
SoBRLPickDetail::getSourceId(void) const
{
    return this->sourceId;
}

void
SoBRLPickDetail::setRegionId(int id)
{
    this->regionId = id;
}

int
SoBRLPickDetail::getRegionId(void) const
{
    return this->regionId;
}

void
SoBRLPickDetail::setAirCode(int air)
{
    this->airCode = air;
}

int
SoBRLPickDetail::getAirCode(void) const
{
    return this->airCode;
}

void
SoBRLPickDetail::setMaterialId(int material)
{
    this->materialId = material;
}

int
SoBRLPickDetail::getMaterialId(void) const
{
    return this->materialId;
}

void
SoBRLPickDetail::setLos(int value)
{
    this->los = value;
}

int
SoBRLPickDetail::getLos(void) const
{
    return this->los;
}

void
SoBRLPickDetail::setMaterialColor(SbBool valid, const SbColor &color)
{
    this->materialColorValid = valid;
    this->materialColor = color;
}

SbBool
SoBRLPickDetail::hasMaterialColor(void) const
{
    return this->materialColorValid;
}

const SbColor &
SoBRLPickDetail::getMaterialColor(void) const
{
    return this->materialColor;
}

void
SoBRLPickDetail::setMaterialShader(const SbString &shader)
{
    this->materialShader = shader;
}

const SbString &
SoBRLPickDetail::getMaterialShader(void) const
{
    return this->materialShader;
}

void
SoBRLPickDetail::setPrimitive(PrimitiveKind kind, int index)
{
    this->primitiveKind = kind;
    this->primitiveIndex = index;
}

SoBRLPickDetail::PrimitiveKind
SoBRLPickDetail::getPrimitiveKind(void) const
{
    return this->primitiveKind;
}

int
SoBRLPickDetail::getPrimitiveIndex(void) const
{
    return this->primitiveIndex;
}

void
SoBRLPickDetail::setEditIntent(const SbString &id, const SbString &role)
{
    this->editIntentId = id;
    this->editIntentRole = role;
}

const SbString &
SoBRLPickDetail::getEditIntentId(void) const
{
    return this->editIntentId;
}

const SbString &
SoBRLPickDetail::getEditIntentRole(void) const
{
    return this->editIntentRole;
}

void
SoBRLPickDetail::setFaceVertexIndices(int indexA, int indexB, int indexC)
{
    this->faceVertexIndex[0] = indexA;
    this->faceVertexIndex[1] = indexB;
    this->faceVertexIndex[2] = indexC;
}

int
SoBRLPickDetail::getFaceVertexIndex(int vertexSlot) const
{
    if (vertexSlot < 0 || vertexSlot >= 3)
	return -1;
    return this->faceVertexIndex[vertexSlot];
}

int
SoBRLPickDetail::getFaceVertexIndexA(void) const
{
    return this->getFaceVertexIndex(0);
}

int
SoBRLPickDetail::getFaceVertexIndexB(void) const
{
    return this->getFaceVertexIndex(1);
}

int
SoBRLPickDetail::getFaceVertexIndexC(void) const
{
    return this->getFaceVertexIndex(2);
}

void
SoBRLPickDetail::setNearestFaceEdge(int edgeSlot, int indexA, int indexB)
{
    this->nearestFaceEdgeSlot = edgeSlot;
    this->nearestFaceEdgeVertexIndex[0] = indexA;
    this->nearestFaceEdgeVertexIndex[1] = indexB;
}

int
SoBRLPickDetail::getNearestFaceEdgeSlot(void) const
{
    return this->nearestFaceEdgeSlot;
}

int
SoBRLPickDetail::getNearestFaceEdgeVertexIndexA(void) const
{
    return this->nearestFaceEdgeVertexIndex[0];
}

int
SoBRLPickDetail::getNearestFaceEdgeVertexIndexB(void) const
{
    return this->nearestFaceEdgeVertexIndex[1];
}

void
SoBRLPickDetail::setNearestFaceVertex(int vertexSlot, int vertexIndex)
{
    this->nearestFaceVertexSlot = vertexSlot;
    this->nearestFaceVertexIndex = vertexIndex;
}

int
SoBRLPickDetail::getNearestFaceVertexSlot(void) const
{
    return this->nearestFaceVertexSlot;
}

int
SoBRLPickDetail::getNearestFaceVertexIndex(void) const
{
    return this->nearestFaceVertexIndex;
}

void
SoBRLPickDetail::setModelPoint(const SbVec3f &point)
{
    this->modelPoint = point;
}

const SbVec3f &
SoBRLPickDetail::getModelPoint(void) const
{
    return this->modelPoint;
}
