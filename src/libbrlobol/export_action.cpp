/*                 E X P O R T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "brlobol/database_source.h"
#include "brlobol/export_action.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/view_lod.h"
#include "brlobol/vlist_shape.h"

#include "bu/path.h"

#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoTransformation.h>

#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_map>
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

static SbString
export_record_display_name(const SbString &path, const SbString &sourceName)
{
    if (sourceName.getLength() > 0)
	return sourceName;
    return path;
}

static uint64_t
export_identity_value(const char *identity)
{
    if (!identity || !identity[0])
	return 0;

    uint64_t hash = 14695981039346656037ULL;
    while (*identity) {
	hash ^= static_cast<unsigned char>(*identity);
	hash *= 1099511628211ULL;
	identity++;
    }

    return hash ? hash : 1;
}

static SbString
export_record_identity_fallback(const SbString &path,
				const SbString &sourceName)
{
    return path.getLength() > 0 ? path : sourceName;
}

template <typename Record>
static void
export_update_record_identity_values(Record &record)
{
    const SbString cacheIdentity = record.cacheIdentity.getLength() > 0 ?
				   record.cacheIdentity : record.sourceIdentity;
    record.cacheIdentityValue =
	export_identity_value(cacheIdentity.getString());
    record.sourceIdentityValue =
	export_identity_value(record.sourceIdentity.getString());
}

template <typename Record>
static void
export_init_record_metadata(Record &record,
			    const SbString &path,
			    const SbString &sourceName,
			    const SbString &sourceType,
			    const char *fallbackGeometryKind)
{
    record.displayName = export_record_display_name(path, sourceName);
    record.geometryName = sourceName.getLength() > 0 ? sourceName : sourceType;
    record.cacheIdentity = "";
    record.sourceIdentity = export_record_identity_fallback(path, sourceName);
    export_update_record_identity_values(record);
    record.databaseIntent = 0;
    record.overlayIntent = 0;
    record.hudIntent = 0;
    record.localSource = 0;
    record.sharedSource = 0;
    record.nonDatabaseSource = 0;
    record.drawMode = 0;
    record.recordRole = "";
    record.geometryKind = fallbackGeometryKind ? fallbackGeometryKind : "";
}

template <typename Record, typename Shape>
static void
export_apply_shape_metadata(Record &record,
			    const Shape *shape,
			    const char *fallbackGeometryKind)
{
    if (!shape)
	return;

    const SbString &displayName = shape->displayName.getValue();
    const SbString &geometryName = shape->geometryName.getValue();
    const SbString &sourceIdentity = shape->sourceIdentity.getValue();
    const SbString &geometryKind = shape->geometryKind.getValue();

    record.displayName = displayName.getLength() > 0 ?
			 displayName : export_record_display_name(shape->sourcePath.getValue(),
			     shape->sourceName.getValue());
    record.geometryName = geometryName.getLength() > 0 ?
			  geometryName : shape->sourceName.getValue();
    record.cacheIdentity = shape->cacheIdentity.getValue();
    record.sourceIdentity = sourceIdentity.getLength() > 0 ?
			    sourceIdentity : export_record_identity_fallback(
				shape->sourcePath.getValue(), shape->sourceName.getValue());
    export_update_record_identity_values(record);
    record.databaseIntent = shape->databaseIntent.getValue() ? 1 : 0;
    record.overlayIntent = shape->overlayIntent.getValue() ? 1 : 0;
    record.hudIntent = shape->hudIntent.getValue() ? 1 : 0;
    record.localSource = shape->localSource.getValue() ? 1 : 0;
    record.sharedSource = shape->sharedSource.getValue() ? 1 : 0;
    record.nonDatabaseSource = shape->nonDatabaseSource.getValue() ? 1 : 0;
    record.drawMode = shape->drawMode.getValue();
    record.recordRole = shape->recordRole.getValue();
    record.geometryKind = geometryKind.getLength() > 0 ?
			  geometryKind : SbString(fallbackGeometryKind ? fallbackGeometryKind : "");
}

template <typename Record>
static void
export_apply_source_request_metadata(Record &record,
				     const BRLObolSourceMeshRequest &request,
				     const char *fallbackGeometryKind)
{
    record.displayName = request.displayName.getLength() > 0 ?
			 request.displayName : export_record_display_name(request.path,
			     request.sourceName);
    record.geometryName = request.geometryName.getLength() > 0 ?
			  request.geometryName : request.sourceName;
    record.cacheIdentity = request.cacheIdentity;
    record.sourceIdentity = request.sourceIdentity.getLength() > 0 ?
			    request.sourceIdentity : export_record_identity_fallback(
				request.path, request.sourceName);
    export_update_record_identity_values(record);
    record.databaseIntent = request.databaseIntent ? 1 : 0;
    record.overlayIntent = request.overlayIntent ? 1 : 0;
    record.hudIntent = request.hudIntent ? 1 : 0;
    record.localSource = request.localSource ? 1 : 0;
    record.sharedSource = request.sharedSource ? 1 : 0;
    record.nonDatabaseSource = request.nonDatabaseSource ? 1 : 0;
    record.drawMode = request.drawMode;
    record.recordRole = request.recordRole;
    record.geometryKind = request.geometryKind.getLength() > 0 ?
			  request.geometryKind : SbString(fallbackGeometryKind ? fallbackGeometryKind : "");
}

static void
export_object_key_append_string(std::string &key, const SbString &value)
{
    const char *s = value.getString();
    size_t len = s ? strlen(s) : 0;
    key += std::to_string(len);
    key += ':';
    if (s && len > 0)
	key.append(s, len);
    key += ';';
}

template <typename Record>
static std::string
export_object_key(const Record &record)
{
    std::string key;
    export_object_key_append_string(key, record.path);
    export_object_key_append_string(key, record.sourceName);
    export_object_key_append_string(key, record.sourceType);
    export_object_key_append_string(key, record.cacheIdentity);
    export_object_key_append_string(key, record.sourceIdentity);
    export_object_key_append_string(key, record.recordRole);
    key += std::to_string(record.sourceId);
    key += ';';
    key += std::to_string(record.drawMode);
    key += ';';
    key += std::to_string(record.databaseIntent ? 1 : 0);
    key += std::to_string(record.overlayIntent ? 1 : 0);
    key += std::to_string(record.hudIntent ? 1 : 0);
    key += std::to_string(record.localSource ? 1 : 0);
    key += std::to_string(record.sharedSource ? 1 : 0);
    key += std::to_string(record.nonDatabaseSource ? 1 : 0);
    return key;
}

template <typename Record>
static void
export_object_init_common(SoBRLExportAction::ObjectRecord &object,
			  const Record &record)
{
    object.sequence = record.sequence;
    object.path = record.path;
    object.sourceName = record.sourceName;
    object.sourceType = record.sourceType;
    object.sourceId = record.sourceId;
    object.displayName = record.displayName;
    object.geometryName = record.geometryName;
    object.cacheIdentity = record.cacheIdentity;
    object.sourceIdentity = record.sourceIdentity;
    object.cacheIdentityValue = record.cacheIdentityValue;
    object.sourceIdentityValue = record.sourceIdentityValue;
    object.databaseIntent = record.databaseIntent;
    object.overlayIntent = record.overlayIntent;
    object.hudIntent = record.hudIntent;
    object.localSource = record.localSource;
    object.sharedSource = record.sharedSource;
    object.nonDatabaseSource = record.nonDatabaseSource;
    object.drawMode = record.drawMode;
    object.recordRole = record.recordRole;
    object.geometryKind = record.geometryKind;
    object.selected = record.selected ? 1 : 0;
    object.highlighted = record.highlighted ? 1 : 0;
    object.visible = 1;
    object.colorOverride = record.colorOverride ? 1 : 0;
    object.color = record.color;
    object.bounds.makeEmpty();
    object.lineIndices.clear();
    object.pointIndices.clear();
    object.triangleIndices.clear();
}

template <typename Record>
static void
export_object_update_common(SoBRLExportAction::ObjectRecord &object,
			    const Record &record)
{
    if (record.sequence < object.sequence)
	object.sequence = record.sequence;
    object.selected = object.selected || record.selected;
    object.highlighted = object.highlighted || record.highlighted;
    object.visible = 1;
    if (!object.colorOverride && record.colorOverride) {
	object.colorOverride = 1;
	object.color = record.color;
    }
    if (object.geometryKind.getLength() == 0) {
	object.geometryKind = record.geometryKind;
    } else if (record.geometryKind.getLength() > 0 &&
	       bu_strcmp(object.geometryKind.getString(),
		      record.geometryKind.getString()) != 0) {
	object.geometryKind = "mixed";
    }
}

static void
export_object_add_line(SoBRLExportAction::ObjectRecord &object,
		       const SoBRLExportAction::LineRecord &record,
		       int index)
{
    export_object_update_common(object, record);
    object.lineIndices.push_back(index);
    object.bounds.extendBy(record.a);
    object.bounds.extendBy(record.b);
}

static void
export_object_add_point(SoBRLExportAction::ObjectRecord &object,
			const SoBRLExportAction::PointRecord &record,
			int index)
{
    export_object_update_common(object, record);
    object.pointIndices.push_back(index);
    object.bounds.extendBy(record.point);
    if (record.pointScaleValid && record.pointScale > 0.0f) {
	const SbVec3f radius(record.pointScale, record.pointScale,
			     record.pointScale);
	object.bounds.extendBy(record.point - radius);
	object.bounds.extendBy(record.point + radius);
    }
}

static void
export_object_add_triangle(SoBRLExportAction::ObjectRecord &object,
			   const SoBRLExportAction::TriangleRecord &record,
			   int index)
{
    export_object_update_common(object, record);
    object.triangleIndices.push_back(index);
    object.bounds.extendBy(record.a);
    object.bounds.extendBy(record.b);
    object.bounds.extendBy(record.c);
}

static int
export_object_record_matches_glob(
    const SoBRLExportAction::ObjectRecord &record,
    const char *glob)
{
    if (!glob || !glob[0])
	return 1;

    return bu_path_match(glob, record.path.getString(), 0) == 0 ||
	   bu_path_match(glob, record.displayName.getString(), 0) == 0 ||
	   bu_path_match(glob, record.sourceName.getString(), 0) == 0;
}

static int
export_object_record_matches_query(
    const SoBRLExportAction::ObjectRecord &record,
    unsigned int queryFlags,
    const char *glob,
    int drawMode)
{
    if ((queryFlags & SoBRLExportAction::QUERY_VISIBLE_ONLY) &&
	!record.visible)
	return 0;
    if ((queryFlags & SoBRLExportAction::QUERY_LOCAL_ONLY) &&
	!record.localSource)
	return 0;
    if (drawMode >= 0 && record.drawMode != drawMode)
	return 0;
    if (!export_object_record_matches_glob(record, glob))
	return 0;

    const int wantsDatabase =
	(queryFlags & SoBRLExportAction::QUERY_DATABASE_OBJECTS) ? 1 : 0;
    const int wantsView =
	(queryFlags & SoBRLExportAction::QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wantsDatabase || wantsView) {
	const int isDatabase = record.databaseIntent &&
			       !record.nonDatabaseSource;
	const int isView = record.overlayIntent || record.hudIntent ||
			   record.localSource || record.sharedSource ||
			   record.nonDatabaseSource;
	if (wantsDatabase && isDatabase)
	    return 1;
	if (wantsView && isView)
	    return 1;
	return 0;
    }

    return 1;
}

static const SoBRLExportAction::TriangleRecord *
export_object_first_triangle(
    const std::vector<SoBRLExportAction::TriangleRecord> &triangles,
    const SoBRLExportAction::ObjectRecord &record)
{
    for (size_t i = 0; i < record.triangleIndices.size(); i++) {
	int triangleIndex = record.triangleIndices[i];
	if (triangleIndex < 0 ||
	    static_cast<size_t>(triangleIndex) >= triangles.size())
	    continue;
	return &triangles[static_cast<size_t>(triangleIndex)];
    }

    return NULL;
}

static int
export_triangle_vertex_id(const SoBRLExportAction::TriangleRecord &triangle,
			  size_t vertexSlot)
{
    if (vertexSlot == 0)
	return triangle.vertexIndexA;
    if (vertexSlot == 1)
	return triangle.vertexIndexB;
    return triangle.vertexIndexC;
}

static SbVec3f
export_triangle_vertex_point(
    const SoBRLExportAction::TriangleRecord &triangle,
    size_t vertexSlot)
{
    if (vertexSlot == 0)
	return triangle.a;
    if (vertexSlot == 1)
	return triangle.b;
    return triangle.c;
}

static int
export_surface_compact_vertex_slot(std::vector<int> &vertexIds,
				   std::vector<SbVec3f> &points, int vertexId, const SbVec3f &point)
{
    if (vertexId >= 0) {
	for (size_t i = 0; i < vertexIds.size(); i++) {
	    if (vertexIds[i] == vertexId)
		return static_cast<int>(i);
	}
    }

    vertexIds.push_back(vertexId);
    points.push_back(point);
    return static_cast<int>(points.size() - 1);
}

static void
export_object_surface_compact_detail(
    const std::vector<SoBRLExportAction::TriangleRecord> &triangles,
    const SoBRLExportAction::ObjectRecord &record,
    std::vector<SbVec3f> *pointsOut,
    std::vector<int> *indicesOut)
{
    std::vector<int> vertexIds;
    std::vector<SbVec3f> points;
    std::vector<int> indices;

    points.reserve(record.triangleIndices.size() * 3);
    indices.reserve(record.triangleIndices.size() * 3);
    for (size_t i = 0; i < record.triangleIndices.size(); i++) {
	int triangleIndex = record.triangleIndices[i];
	if (triangleIndex < 0 ||
	    static_cast<size_t>(triangleIndex) >= triangles.size())
	    continue;
	const SoBRLExportAction::TriangleRecord &triangle =
	    triangles[static_cast<size_t>(triangleIndex)];
	for (size_t vertexSlot = 0; vertexSlot < 3; vertexSlot++) {
	    indices.push_back(export_surface_compact_vertex_slot(vertexIds,
			      points,
			      export_triangle_vertex_id(triangle, vertexSlot),
			      export_triangle_vertex_point(triangle, vertexSlot)));
	}
    }

    if (pointsOut)
	pointsOut->swap(points);
    if (indicesOut)
	indicesOut->swap(indices);
}

SoBRLExportAction::SoBRLExportAction(void) :
    lineCount(0),
    pointCount(0),
    triangleCount(0),
    geometryPolicy(FULL_DETAIL),
    recordStorageEnabled(TRUE),
    skippedLodDisplayMeshCount(0),
    recordSequence(0),
    objectRecordsDirty(TRUE)
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
    SO_ENABLE(SoBRLExportAction, SoBRLViewLodElement);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLExportAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLDatabaseSource,
			 SoBRLExportAction::databaseSourceAction);
    SO_ACTION_ADD_METHOD(SoBRLVListShape, SoBRLExportAction::vlistShapeAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape, SoBRLExportAction::meshShapeAction);
}

uint64_t
SoBRLExportAction::identityValue(const char *identity)
{
    return export_identity_value(identity);
}

uint64_t
SoBRLExportAction::identityValue(const SbString &identity)
{
    return export_identity_value(identity.getString());
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

int
SoBRLExportAction::getObjectRecordCount(void) const
{
    this->rebuildObjectRecords();
    return static_cast<int>(this->objects.size());
}

const SoBRLExportAction::ObjectRecord &
SoBRLExportAction::getObjectRecord(int index) const
{
    this->rebuildObjectRecords();
    return this->objects.at(static_cast<size_t>(index));
}

int
SoBRLExportAction::collectObjectRecords(std::vector<ObjectRecord> &records,
					unsigned int queryFlags,
					const char *glob,
					int drawMode) const
{
    this->rebuildObjectRecords();
    records.clear();
    for (size_t i = 0; i < this->objects.size(); i++) {
	if (!export_object_record_matches_query(this->objects[i],
						queryFlags, glob, drawMode))
	    continue;
	records.push_back(this->objects[i]);
    }
    return static_cast<int>(records.size());
}

SbBool
SoBRLExportAction::getObjectRecordLineSummary(const ObjectRecord &record,
	ObjectLineSummary &summary) const
{
    summary.valid = 0;
    summary.pointCount = 0;
    summary.segmentCount = 0;
    summary.cacheIdentity = "";
    summary.sourceIdentity = "";
    summary.cacheIdentityValue = 0;
    summary.sourceIdentityValue = 0;
    if (record.lineIndices.empty())
	return FALSE;

    summary.valid = 1;
    summary.segmentCount = record.lineIndices.size();
    summary.pointCount = summary.segmentCount * 2;
    summary.cacheIdentity = record.cacheIdentity;
    summary.sourceIdentity = record.sourceIdentity;
    summary.cacheIdentityValue = record.cacheIdentityValue;
    summary.sourceIdentityValue = record.sourceIdentityValue;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordLinePoint(const ObjectRecord &record,
	size_t index,
	SbVec3f &point) const
{
    point = SbVec3f(0.0f, 0.0f, 0.0f);
    size_t segmentIndex = index / 2;
    if (segmentIndex >= record.lineIndices.size())
	return FALSE;
    int lineIndex = record.lineIndices[segmentIndex];
    if (lineIndex < 0 ||
	static_cast<size_t>(lineIndex) >= this->lines.size())
	return FALSE;

    const LineRecord &line = this->lines[static_cast<size_t>(lineIndex)];
    point = (index % 2) ? line.b : line.a;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordLineCommand(const ObjectRecord &record,
	size_t index,
	int &command) const
{
    size_t segmentIndex = index / 2;
    if (segmentIndex >= record.lineIndices.size())
	return FALSE;
    int lineIndex = record.lineIndices[segmentIndex];
    if (lineIndex < 0 ||
	static_cast<size_t>(lineIndex) >= this->lines.size())
	return FALSE;
    command = (index % 2) ? LINE_DRAW : LINE_MOVE;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordPointSummary(const ObjectRecord &record,
	ObjectPointSummary &summary) const
{
    summary.valid = 0;
    summary.pointCount = 0;
    summary.cacheIdentity = "";
    summary.sourceIdentity = "";
    summary.cacheIdentityValue = 0;
    summary.sourceIdentityValue = 0;
    if (record.pointIndices.empty())
	return FALSE;

    summary.valid = 1;
    summary.pointCount = record.pointIndices.size();
    summary.cacheIdentity = record.cacheIdentity;
    summary.sourceIdentity = record.sourceIdentity;
    summary.cacheIdentityValue = record.cacheIdentityValue;
    summary.sourceIdentityValue = record.sourceIdentityValue;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordPoint(const ObjectRecord &record,
					size_t index,
					SbVec3f &point) const
{
    point = SbVec3f(0.0f, 0.0f, 0.0f);
    if (index >= record.pointIndices.size())
	return FALSE;
    int pointIndex = record.pointIndices[index];
    if (pointIndex < 0 ||
	static_cast<size_t>(pointIndex) >= this->points.size())
	return FALSE;

    point = this->points[static_cast<size_t>(pointIndex)].point;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordSurfaceSummary(const ObjectRecord &record,
	ObjectSurfaceSummary &summary) const
{
    summary.valid = 0;
    summary.pointCount = 0;
    summary.normalCount = 0;
    summary.indexCount = 0;
    summary.faceCount = 0;
    summary.normalsPerIndex = 0;
    summary.materialColorValid = 0;
    summary.materialColor = SbColor(0.0f, 0.0f, 0.0f);
    summary.materialDrawMode = record.drawMode;
    summary.materialHighlighted = record.highlighted;
    summary.cacheIdentity = "";
    summary.sourceIdentity = "";
    summary.cacheIdentityValue = 0;
    summary.sourceIdentityValue = 0;
    std::vector<SbVec3f> surfacePoints;
    std::vector<int> surfaceIndices;
    export_object_surface_compact_detail(this->triangles, record,
					 &surfacePoints, &surfaceIndices);
    const TriangleRecord *firstTriangle =
	export_object_first_triangle(this->triangles, record);
    if (!firstTriangle)
	return FALSE;

    summary.valid = 1;
    summary.faceCount = surfaceIndices.size() / 3;
    summary.indexCount = surfaceIndices.size();
    summary.pointCount = surfacePoints.size();
    summary.normalCount = firstTriangle->lodNormalCount;
    summary.normalsPerIndex = firstTriangle->lodHasNormals ? 1 : 0;
    summary.materialColorValid = firstTriangle->materialColorValid;
    summary.materialColor = firstTriangle->materialColor;
    summary.materialDrawMode = firstTriangle->drawMode;
    summary.materialHighlighted = record.highlighted;
    summary.cacheIdentity = record.cacheIdentity;
    summary.sourceIdentity = record.sourceIdentity;
    summary.cacheIdentityValue = record.cacheIdentityValue;
    summary.sourceIdentityValue = record.sourceIdentityValue;
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordSurfacePoint(const ObjectRecord &record,
	size_t index,
	SbVec3f &point) const
{
    point = SbVec3f(0.0f, 0.0f, 0.0f);

    std::vector<SbVec3f> surfacePoints;
    export_object_surface_compact_detail(this->triangles, record,
					 &surfacePoints, NULL);
    if (index >= surfacePoints.size())
	return FALSE;

    point = surfacePoints[index];
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordSurfaceDetail(const ObjectRecord &record,
	std::vector<SbVec3f> *surfacePoints,
	std::vector<int> *surfaceIndices) const
{
    if (record.triangleIndices.empty()) {
	if (surfacePoints)
	    surfacePoints->clear();
	if (surfaceIndices)
	    surfaceIndices->clear();
	return FALSE;
    }

    export_object_surface_compact_detail(this->triangles, record,
					 surfacePoints, surfaceIndices);
    return TRUE;
}

SbBool
SoBRLExportAction::getObjectRecordSurfaceIndex(const ObjectRecord &record,
	size_t index,
	int &out) const
{
    std::vector<int> surfaceIndices;
    export_object_surface_compact_detail(this->triangles, record, NULL,
					 &surfaceIndices);
    if (index >= surfaceIndices.size())
	return FALSE;

    out = surfaceIndices[index];
    return TRUE;
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
	export_apply_source_request_metadata(this->triangles.back(),
					     sourceRequest, "surface");
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
    if (node->isOfType(SoGroup::getClassTypeId()) ||
	node->isOfType(SoTransformation::getClassTypeId()))
	node->doAction(action);
}

void
SoBRLExportAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLExportAction *exportAction =
	static_cast<SoBRLExportAction *>(action);
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
    const SbMatrix &parentToWorld =
	SoModelMatrixElement::get(action->getState());
    if (source->exportCompactInstances(exportAction, parentToWorld))
	return;

    source->doAction(action);
}

void
SoBRLExportAction::vlistShapeAction(SoAction *action, SoNode *node)
{
    SoBRLExportAction *exportAction = static_cast<SoBRLExportAction *>(action);
    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    if (!shape->visible.getValue())
	return;

    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (exportAction->recordStorageEnabled) {
	export_reserve_records(exportAction->lines,
			       static_cast<size_t>(shape->getSegmentCount()));
	export_reserve_records(exportAction->points,
			       static_cast<size_t>(shape->getPointPrimitiveCount()));
    }

    SbVec3f last;
    SbBool haveLast = FALSE;
    int segmentIndex = 0;
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
	    if (haveLast) {
		SbVec3f worldA;
		SbVec3f worldB;
		localToWorld.multVecMatrix(last, worldA);
		localToWorld.multVecMatrix(geom->point[i], worldB);

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
		    export_apply_shape_metadata(exportAction->lines.back(),
						shape, "line");
		}
		segmentIndex++;
	    }
	    last = geom->point[i];
	    haveLast = TRUE;
	    continue;
	}

	if (geom->command[i] == SoBRLVListShape::POINT) {
	    SbVec3f worldPoint;
	    localToWorld.multVecMatrix(geom->point[i], worldPoint);

	    float pointScale = 0.0f;
	    const int pointScaleValid =
		shape->getPointScale(i, pointScale) ? 1 : 0;
	    if (pointScaleValid)
		pointScale = export_transform_point_scale(localToWorld,
			     geom->point[i], worldPoint, pointScale);

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
	    export_apply_shape_metadata(exportAction->points.back(),
					shape, "point");
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

    const BRLObolViewLodState::MeshPayload *viewPayload =
	brlobol_view_lod_mesh_for_action(action, shape);
    const SbBool useFullDetail =
	(exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	 shape->hasFullDetailMesh()) ? TRUE : FALSE;
    const SbBool useSourceBackedFullDetail =
	(exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	 shape->needsSourceBackedFullDetail()) ? TRUE : FALSE;
    if (exportAction->geometryPolicy == SoBRLExportAction::FULL_DETAIL &&
	(shape->isLodDisplayActive() || viewPayload))
	exportAction->skippedLodDisplayMeshCount++;
    const SbMatrix &localToWorld = SoModelMatrixElement::get(action->getState());
    if (useSourceBackedFullDetail) {
	exportAction->appendSourceBackedFullDetailRequest(shape, localToWorld);
	return;
    }

    const SbBool useViewPayload =
	(!useFullDetail &&
	 exportAction->geometryPolicy == SoBRLExportAction::DISPLAY_LEVEL &&
	 viewPayload) ? TRUE : FALSE;
    int triangleCount = useFullDetail ?
			shape->getFullDetailTriangleCount() :
			(useViewPayload ? viewPayload->getTriangleCount() :
			 shape->getTriangleCount());
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
	} else if (useViewPayload) {
	    if (!viewPayload->getTriangle(i, a, b, c))
		continue;
	    if (!viewPayload->getTriangleVertexIndices(i, ia, ib, ic))
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

	const SbBool lodAvailable = useViewPayload ? TRUE :
				    shape->lodAvailable.getValue();
	const int lodActiveLevel = useViewPayload ? viewPayload->activeLevel :
				   shape->lodActiveLevel.getValue();
	const uint32_t lodFaceCount = useViewPayload ?
				      static_cast<uint32_t>(viewPayload->counts.faceCount > UINT32_MAX ?
					  UINT32_MAX : viewPayload->counts.faceCount) :
				      shape->lodFaceCount.getValue();
	const uint32_t lodPointCount = useViewPayload ?
				       static_cast<uint32_t>(viewPayload->counts.pointCount > UINT32_MAX ?
					   UINT32_MAX : viewPayload->counts.pointCount) :
				       shape->lodPointCount.getValue();
	const uint32_t lodOriginalPointCount = useViewPayload ?
					       static_cast<uint32_t>(
						       viewPayload->counts.originalPointCount > UINT32_MAX ?
						       UINT32_MAX : viewPayload->counts.originalPointCount) :
					       shape->lodOriginalPointCount.getValue();
	const uint32_t lodNormalCount = useViewPayload ?
					static_cast<uint32_t>(viewPayload->counts.normalCount > UINT32_MAX ?
					    UINT32_MAX : viewPayload->counts.normalCount) :
					shape->lodNormalCount.getValue();
	const SbVec3f lodBoundsMin =
	    useViewPayload && !viewPayload->bounds.isEmpty() ?
	    viewPayload->bounds.getMin() : shape->lodBoundsMin.getValue();
	const SbVec3f lodBoundsMax =
	    useViewPayload && !viewPayload->bounds.isEmpty() ?
	    viewPayload->bounds.getMax() : shape->lodBoundsMax.getValue();

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
				     lodAvailable,
				     lodActiveLevel,
				     lodFaceCount,
				     lodPointCount,
				     lodOriginalPointCount,
				     lodNormalCount,
				     useViewPayload ? viewPayload->hasSnappedPoints :
				     shape->lodHasSnappedPoints.getValue(),
				     useViewPayload ? viewPayload->hasNormals :
				     shape->lodHasNormals.getValue(),
				     lodBoundsMin,
				     lodBoundsMax,
				     shape->colorOverride.getValue(),
				     shape->color.getValue(), worldA, worldB, worldC);
	export_apply_shape_metadata(exportAction->triangles.back(),
				    shape, "surface");
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
    this->recordSequence = 0;
    this->objects.clear();
    this->objectRecordsDirty = TRUE;
}

void
SoBRLExportAction::invalidateObjectRecords(void)
{
    this->objectRecordsDirty = TRUE;
}

struct ExportPrimitiveRef {
    size_t sequence;
    int kind;
    int index;
};

void
SoBRLExportAction::rebuildObjectRecords(void) const
{
    if (!this->objectRecordsDirty)
	return;

    this->objects.clear();
    std::vector<ExportPrimitiveRef> refs;
    refs.reserve(this->lines.size() + this->points.size() +
		 this->triangles.size());
    for (size_t i = 0; i < this->lines.size(); i++) {
	ExportPrimitiveRef ref = {this->lines[i].sequence, 0,
				  static_cast<int>(i)
				 };
	refs.push_back(ref);
    }
    for (size_t i = 0; i < this->points.size(); i++) {
	ExportPrimitiveRef ref = {this->points[i].sequence, 1,
				  static_cast<int>(i)
				 };
	refs.push_back(ref);
    }
    for (size_t i = 0; i < this->triangles.size(); i++) {
	ExportPrimitiveRef ref = {this->triangles[i].sequence, 2,
				  static_cast<int>(i)
				 };
	refs.push_back(ref);
    }
    std::sort(refs.begin(), refs.end(),
    [](const ExportPrimitiveRef &a, const ExportPrimitiveRef &b) {
	return a.sequence < b.sequence;
    });

    std::unordered_map<std::string, size_t> objectIndexByKey;
    objectIndexByKey.reserve(refs.size());
    for (size_t i = 0; i < refs.size(); i++) {
	const ExportPrimitiveRef &ref = refs[i];
	std::string key;
	if (ref.kind == 0)
	    key = export_object_key(this->lines[static_cast<size_t>(ref.index)]);
	else if (ref.kind == 1)
	    key = export_object_key(this->points[static_cast<size_t>(ref.index)]);
	else
	    key = export_object_key(this->triangles[static_cast<size_t>(ref.index)]);

	size_t objectIndex = 0;
	std::unordered_map<std::string, size_t>::const_iterator found =
	    objectIndexByKey.find(key);
	if (found == objectIndexByKey.end()) {
	    ObjectRecord object;
	    if (ref.kind == 0)
		export_object_init_common(object,
					  this->lines[static_cast<size_t>(ref.index)]);
	    else if (ref.kind == 1)
		export_object_init_common(object,
					  this->points[static_cast<size_t>(ref.index)]);
	    else
		export_object_init_common(object,
					  this->triangles[static_cast<size_t>(ref.index)]);
	    this->objects.push_back(object);
	    objectIndex = this->objects.size() - 1;
	    objectIndexByKey[key] = objectIndex;
	} else {
	    objectIndex = found->second;
	}

	ObjectRecord &object = this->objects[objectIndex];
	if (ref.kind == 0)
	    export_object_add_line(object,
				   this->lines[static_cast<size_t>(ref.index)], ref.index);
	else if (ref.kind == 1)
	    export_object_add_point(object,
				    this->points[static_cast<size_t>(ref.index)], ref.index);
	else
	    export_object_add_triangle(object,
				       this->triangles[static_cast<size_t>(ref.index)], ref.index);
    }

    std::sort(this->objects.begin(), this->objects.end(),
    [](const ObjectRecord &a, const ObjectRecord &b) {
	return a.sequence < b.sequence;
    });
    this->objectRecordsDirty = FALSE;
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
    record.sequence = this->recordSequence++;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    export_init_record_metadata(record, path, sourceName, sourceType, "line");
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
    this->invalidateObjectRecords();
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
    record.sequence = this->recordSequence++;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    export_init_record_metadata(record, path, sourceName, sourceType, "point");
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
    this->invalidateObjectRecords();
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
    record.sequence = this->recordSequence++;
    record.path = path;
    record.sourceName = sourceName;
    record.sourceType = sourceType;
    record.sourceId = sourceId;
    export_init_record_metadata(record, path, sourceName, sourceType, "surface");
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
    this->invalidateObjectRecords();
    this->bounds.extendBy(a);
    this->bounds.extendBy(b);
    this->bounds.extendBy(c);
}
