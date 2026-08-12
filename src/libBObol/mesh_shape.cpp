/*                   M E S H _ S H A P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLodRealization.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BViewLod.h"

#include <Inventor/SbVec3f.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoFaceDetail.h>
#include <Inventor/details/SoPointDetail.h>
#include <Inventor/elements/SoGLCacheContextElement.h>
#include <Inventor/elements/SoGLDisplayList.h>
#include <Inventor/gl.h>

#include <limits>
#include <map>
#include <sstream>
#include <stdint.h>
#include <string>

SO_NODE_SOURCE(SoBRLMeshShape);

static uint32_t
mesh_lod_count_to_field(uint64_t count)
{
    return count > UINT32_MAX ? UINT32_MAX : static_cast<uint32_t>(count);
}

static SbString
mesh_lod_uint64_string(uint64_t value)
{
    std::ostringstream out;
    out << value;
    return SbString(out.str().c_str());
}

SoBRLMeshShape::SoBRLMeshShape(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLMeshShape);

    this->sourceMeshFaceCount = 0;
    this->sourceMeshPointCount = 0;
    this->sourceMeshBounds.makeEmpty();
    this->sourceMeshMetricsValid = FALSE;

    SO_NODE_ADD_EMPTY_MFIELD(point);
    SO_NODE_ADD_EMPTY_MFIELD(coordIndex);
    SO_NODE_ADD_EMPTY_MFIELD(normal);
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
    SO_NODE_ADD_FIELD(lodAvailable, (FALSE));
    SO_NODE_ADD_FIELD(lodActiveCut, (-1));
    SO_NODE_ADD_FIELD(lodFaceCount, (0));
    SO_NODE_ADD_FIELD(lodPointCount, (0));
    SO_NODE_ADD_FIELD(lodOriginalPointCount, (0));
    SO_NODE_ADD_FIELD(lodNormalCount, (0));
    SO_NODE_ADD_FIELD(lodHasSnappedPoints, (FALSE));
    SO_NODE_ADD_FIELD(lodHasNormals, (FALSE));
    SO_NODE_ADD_FIELD(lodBoundsMin, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodBoundsMax, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodStagedAvailable, (FALSE));
    SO_NODE_ADD_FIELD(lodResultKind, (BOBOL_LOD_RESULT_NONE));
    SO_NODE_ADD_FIELD(lodQualityTier, (BOBOL_LOD_QUALITY_METADATA));
    SO_NODE_ADD_FIELD(lodProviderStatus, (BOBOL_LOD_PROVIDER_UNKNOWN));
    SO_NODE_ADD_FIELD(lodCacheKey, (""));
    SO_NODE_ADD_FIELD(lodProviderId, (""));
    SO_NODE_ADD_FIELD(lodProviderVersion, (""));
    SO_NODE_ADD_FIELD(lodDiagnostic, (""));
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencyPath);
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencyName);
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencySourceRevision);
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencySourceHash);
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencyQualityTier);
    SO_NODE_ADD_EMPTY_MFIELD(lodDependencyOptional);
    SO_NODE_ADD_EMPTY_MFIELD(lodAttributeName);
    SO_NODE_ADD_EMPTY_MFIELD(lodAttributeValue);
    SO_NODE_ADD_FIELD(lodProxyKind, (BOBOL_LOD_PROXY_NONE));
    SO_NODE_ADD_FIELD(lodProxyCenter, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodProxyAxisX, (SbVec3f(1.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodProxyAxisY, (SbVec3f(0.0f, 1.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodProxyAxisZ, (SbVec3f(0.0f, 0.0f, 1.0f)));
    SO_NODE_ADD_FIELD(lodProxyHalfExtents, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodProxyRadius, (0.0f));
    SO_NODE_ADD_FIELD(lodBacked, (FALSE));
    SO_NODE_ADD_FIELD(lodDisplayActive, (FALSE));
    SO_NODE_ADD_FIELD(lodPreserveFullDetail, (TRUE));
    SO_NODE_ADD_FIELD(pickGeometryPolicy, (PICK_DISPLAY_LEVEL));
    SO_NODE_ADD_EMPTY_MFIELD(selectedPrimitive);
    SO_NODE_ADD_EMPTY_MFIELD(highlightedPrimitive);
    SO_NODE_ADD_FIELD(sharedGeometry, (NULL));
}

SoBRLMeshShape::~SoBRLMeshShape(void)
{
    this->clearRenderLists();
}

void
SoBRLMeshShape::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLMeshShape, SoShape, "Shape");
}

void
SoBRLMeshShape::clearRenderLists(void)
{
    for (std::map<int, SoGLDisplayList *>::iterator it =
	 this->renderLists.begin(); it != this->renderLists.end(); ++it) {
	if (it->second)
	    it->second->unref();
    }
    this->renderLists.clear();
    this->renderListSignature.clear();
}

void
SoBRLMeshShape::notify(SoNotList *list)
{
    /* Any field change (new points/indices/normals, mode, selection) may
     * invalidate baked geometry.  Discard the cached display lists; they are
     * rebuilt lazily on the next render.  The view-local LoD payload path does
     * not touch node fields, so it is guarded separately by the render-list
     * signature in mesh_shape_call_cached_geometry(). */
    this->clearRenderLists();
    inherited::notify(list);
}

void
SoBRLMeshShape::setSharedGeometry(SoBRLMeshShape *shape)
{
    this->sharedGeometry = shape;
}

SoBRLMeshShape *
SoBRLMeshShape::getSharedGeometrySource(void)
{
    SoNode *node = this->sharedGeometry.getValue();
    if (node && node != this &&
	node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<SoBRLMeshShape *>(node);
    return this;
}

const SoBRLMeshShape *
SoBRLMeshShape::getSharedGeometrySource(void) const
{
    const SoNode *node = this->sharedGeometry.getValue();
    if (node && node != this &&
	node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node);
    return this;
}

SoBRLMeshShape *
SoBRLMeshShape::getGeometrySource(void)
{
    if (this->getSharedGeometrySource() != this &&
	this->lodDisplayActive.getValue() &&
	this->point.getNum() > 0 &&
	this->coordIndex.getNum() > 0)
	return this;

    return this->getSharedGeometrySource();
}

const SoBRLMeshShape *
SoBRLMeshShape::getGeometrySource(void) const
{
    if (this->getSharedGeometrySource() != this &&
	this->lodDisplayActive.getValue() &&
	this->point.getNum() > 0 &&
	this->coordIndex.getNum() > 0)
	return this;

    return this->getSharedGeometrySource();
}

void
SoBRLMeshShape::setIndexedTriangles(const SbVec3f *points, int pointCount,
				    const int32_t *indices, int indexCount)
{
    this->sharedGeometry = NULL;
    this->clearFullDetailMesh();
    this->lodDisplayActive = FALSE;
    this->setIndexedTriangleFields(points, pointCount, indices, indexCount);
    this->updateSourceMeshMetricsFromFields();
}

void
SoBRLMeshShape::setIndexedTriangles(const SbVec3f *points, int pointCount,
				    const int32_t *indices, int indexCount,
				    const SbVec3f *normals, int normalCount)
{
    this->sharedGeometry = NULL;
    this->clearFullDetailMesh();
    this->lodDisplayActive = FALSE;
    this->setIndexedTriangleFields(points, pointCount, indices, indexCount,
				   normals, normalCount);
    this->updateSourceMeshMetricsFromFields();
}

void
SoBRLMeshShape::setIndexedTriangleFields(const SbVec3f *points, int pointCount,
	const int32_t *indices, int indexCount,
	const SbVec3f *normals, int normalCount)
{
    this->point.setNum(0);
    this->coordIndex.setNum(0);
    this->normal.setNum(0);
    if (!points || pointCount <= 0 || !indices || indexCount <= 0)
	return;

    this->point.setValues(0, pointCount, points);
    this->coordIndex.setValues(0, indexCount, indices);
    if (normals && normalCount == indexCount)
	this->normal.setValues(0, normalCount, normals);
}

int
SoBRLMeshShape::getTriangleCount(void) const
{
    const SoBRLMeshShape *geom = this->getGeometrySource();
    return geom->coordIndex.getNum() / 3;
}

SbBool
SoBRLMeshShape::getTriangle(int triangleIndex, SbVec3f &a, SbVec3f &b, SbVec3f &c) const
{
    int ia = -1;
    int ib = -1;
    int ic = -1;
    if (!this->getTriangleVertexIndices(triangleIndex, ia, ib, ic))
	return FALSE;

    const SoBRLMeshShape *geom = this->getGeometrySource();
    a = geom->point[ia];
    b = geom->point[ib];
    c = geom->point[ic];
    return TRUE;
}

SbBool
SoBRLMeshShape::getTriangleVertexIndices(int triangleIndex, int &indexA, int &indexB, int &indexC) const
{
    indexA = -1;
    indexB = -1;
    indexC = -1;
    if (triangleIndex < 0)
	return FALSE;

    const SoBRLMeshShape *geom = this->getGeometrySource();
    int base = triangleIndex * 3;
    if (base + 2 >= geom->coordIndex.getNum())
	return FALSE;

    indexA = geom->coordIndex[base];
    indexB = geom->coordIndex[base + 1];
    indexC = geom->coordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	indexA >= geom->point.getNum() ||
	indexB >= geom->point.getNum() ||
	indexC >= geom->point.getNum()) {
	indexA = -1;
	indexB = -1;
	indexC = -1;
	return FALSE;
    }

    return TRUE;
}

SbBool
SoBRLMeshShape::hasFullDetailMesh(void) const
{
    const SoBRLMeshShape *geom = this->getGeometrySource();
    return (!geom->fullDetailPoint.empty() &&
	    !geom->fullDetailCoordIndex.empty()) ? TRUE : FALSE;
}

int
SoBRLMeshShape::getFullDetailTriangleCount(void) const
{
    const SoBRLMeshShape *geom = this->getGeometrySource();
    return static_cast<int>(geom->fullDetailCoordIndex.size() / 3);
}

SbBool
SoBRLMeshShape::getFullDetailTriangle(int triangleIndex,
				      SbVec3f &a,
				      SbVec3f &b,
				      SbVec3f &c) const
{
    int ia = -1;
    int ib = -1;
    int ic = -1;
    if (!this->getFullDetailTriangleVertexIndices(triangleIndex, ia, ib, ic))
	return FALSE;

    const SoBRLMeshShape *geom = this->getGeometrySource();
    a = geom->fullDetailPoint[static_cast<size_t>(ia)];
    b = geom->fullDetailPoint[static_cast<size_t>(ib)];
    c = geom->fullDetailPoint[static_cast<size_t>(ic)];
    return TRUE;
}

SbBool
SoBRLMeshShape::getFullDetailTriangleVertexIndices(int triangleIndex,
	int &indexA,
	int &indexB,
	int &indexC) const
{
    indexA = -1;
    indexB = -1;
    indexC = -1;
    if (triangleIndex < 0)
	return FALSE;

    const SoBRLMeshShape *geom = this->getGeometrySource();
    size_t base = static_cast<size_t>(triangleIndex) * 3;
    if (base + 2 >= geom->fullDetailCoordIndex.size())
	return FALSE;

    indexA = geom->fullDetailCoordIndex[base];
    indexB = geom->fullDetailCoordIndex[base + 1];
    indexC = geom->fullDetailCoordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	static_cast<size_t>(indexA) >= geom->fullDetailPoint.size() ||
	static_cast<size_t>(indexB) >= geom->fullDetailPoint.size() ||
	static_cast<size_t>(indexC) >= geom->fullDetailPoint.size()) {
	indexA = -1;
	indexB = -1;
	indexC = -1;
	return FALSE;
    }

    return TRUE;
}

SbBool
SoBRLMeshShape::makeSourceMeshRequest(BObolSourceMeshRequest &request) const
{
    request.clear();
    request.path = this->sourcePath.getValue();
    request.sourceName = this->sourceName.getValue();
    request.sourceType = this->sourceType.getValue();
    request.sourceId = this->sourceId.getValue();
    request.displayName = this->displayName.getValue();
    request.geometryName = this->geometryName.getValue();
    request.cacheIdentity = this->cacheIdentity.getValue();
    request.sourceIdentity = this->sourceIdentity.getValue();
    request.ownerSourceInstanceKey = this->ownerSourceInstanceKey.getValue();
    request.databaseIntent = this->databaseIntent.getValue();
    request.overlayIntent = this->overlayIntent.getValue();
    request.hudIntent = this->hudIntent.getValue();
    request.localSource = this->localSource.getValue();
    request.sharedSource = this->sharedSource.getValue();
    request.nonDatabaseSource = this->nonDatabaseSource.getValue();
    request.drawMode = this->drawMode.getValue();
    request.recordRole = this->recordRole.getValue();
    request.geometryKind = this->geometryKind.getValue();
    request.regionId = this->regionId.getValue();
    request.airCode = this->airCode.getValue();
    request.materialId = this->materialId.getValue();
    request.los = this->los.getValue();
    request.materialColorValid = this->materialColorValid.getValue();
    request.materialColor = this->materialColor.getValue();
    request.materialShader = this->materialShader.getValue();
    request.selected = this->selected.getValue();
    request.highlighted = this->highlighted.getValue();
    request.ghosted = this->ghosted.getValue();
    request.hiddenLine = this->hiddenLine.getValue();
    request.editEmphasis = this->editEmphasis.getValue();
    request.editIntentId = this->editIntentId.getValue();
    request.editIntentRole = this->editIntentRole.getValue();
    request.lodPolicy = this->lodPolicy.getValue();
    request.lodAvailable = this->lodAvailable.getValue();
    request.lodActiveCut = this->lodActiveCut.getValue();
    request.lodFaceCount = this->lodFaceCount.getValue();
    request.lodPointCount = this->lodPointCount.getValue();
    request.lodOriginalPointCount = this->lodOriginalPointCount.getValue();
    request.lodNormalCount = this->lodNormalCount.getValue();
    request.lodHasSnappedPoints = this->lodHasSnappedPoints.getValue();
    request.lodHasNormals = this->lodHasNormals.getValue();
    request.lodBoundsMin = this->lodBoundsMin.getValue();
    request.lodBoundsMax = this->lodBoundsMax.getValue();
    request.colorOverride = this->colorOverride.getValue();
    request.color = this->color.getValue();
    request.transparency = this->transparency.getValue();

    const SoBRLMeshShape *geom = this->getSharedGeometrySource();
    if (!geom->fullDetailPoint.empty() &&
	!geom->fullDetailCoordIndex.empty()) {
	request.faceCount =
	    static_cast<uint64_t>(geom->fullDetailCoordIndex.size() / 3);
	request.pointCount = static_cast<uint64_t>(geom->fullDetailPoint.size());
	for (size_t i = 0; i < geom->fullDetailPoint.size(); i++)
	    request.bounds.extendBy(geom->fullDetailPoint[i]);
    } else if (geom->sourceMeshMetricsValid) {
	request.faceCount = geom->sourceMeshFaceCount;
	request.pointCount = geom->sourceMeshPointCount;
	request.bounds = geom->sourceMeshBounds;
    } else {
	request.faceCount = static_cast<uint64_t>(geom->getTriangleCount());
	request.pointCount = static_cast<uint64_t>(geom->point.getNum());
	for (int i = 0; i < geom->point.getNum(); i++)
	    request.bounds.extendBy(geom->point[i]);
    }

    return (request.faceCount > 0 &&
	    request.pointCount > 0 &&
	    !request.bounds.isEmpty()) ? TRUE : FALSE;
}

SbBool
SoBRLMeshShape::needsSourceBackedFullDetail(void) const
{
    if (this->hasFullDetailMesh())
	return FALSE;

    if (!this->lodDisplayActive.getValue() &&
	!this->lodBacked.getValue() &&
	!(this->lodAvailable.getValue() &&
	  this->lodResultKind.getValue() == BOBOL_LOD_RESULT_MESH))
	return FALSE;

    BObolSourceMeshRequest request;
    return this->makeSourceMeshRequest(request);
}

size_t
SoBRLMeshShape::estimateDisplayMeshBytes(void) const
{
    const SoBRLMeshShape *geom = this->getGeometrySource();
    size_t bytes = static_cast<size_t>(geom->point.getNum()) * sizeof(SbVec3f);
    bytes += static_cast<size_t>(geom->coordIndex.getNum()) * sizeof(int32_t);
    bytes += static_cast<size_t>(geom->normal.getNum()) * sizeof(SbVec3f);
    return bytes;
}

size_t
SoBRLMeshShape::estimateFullDetailMeshBytes(void) const
{
    const SoBRLMeshShape *geom = this->getGeometrySource();
    size_t bytes = geom->fullDetailPoint.capacity() * sizeof(SbVec3f);
    bytes += geom->fullDetailCoordIndex.capacity() * sizeof(int32_t);
    return bytes;
}

size_t
SoBRLMeshShape::estimateResidentMeshBytes(void) const
{
    return this->estimateDisplayMeshBytes() +
	   this->estimateFullDetailMeshBytes();
}

size_t
SoBRLMeshShape::evictFullDetailMesh(void)
{
    if (this->getGeometrySource() != this)
	return 0;

    size_t freedBytes = this->estimateFullDetailMeshBytes();
    if (freedBytes == 0)
	return 0;

    if (this->hasFullDetailMesh())
	this->updateSourceMeshMetricsFromFullDetail();
    this->clearFullDetailMesh();
    return freedBytes;
}

size_t
SoBRLMeshShape::evictDisplayMeshPreservingSourceMetrics(void)
{
    if (this->getGeometrySource() != this)
	return 0;

    size_t freedBytes = this->estimateDisplayMeshBytes();
    if (freedBytes == 0)
	return 0;

    if (!this->sourceMeshMetricsValid)
	this->updateSourceMeshMetricsFromFields();

    this->setIndexedTriangleFields(NULL, 0, NULL, 0);
    this->lodDisplayActive = FALSE;
    return freedBytes;
}

size_t
SoBRLMeshShape::evictActiveDisplayMesh(void)
{
    if (this->getGeometrySource() != this)
	return 0;

    if (!this->lodDisplayActive.getValue())
	return 0;

    size_t freedBytes = this->estimateDisplayMeshBytes();
    if (this->hasFullDetailMesh())
	this->updateSourceMeshMetricsFromFullDetail();
    else if (!this->sourceMeshMetricsValid)
	this->updateSourceMeshMetricsFromFields();

    this->setIndexedTriangleFields(NULL, 0, NULL, 0);
    this->lodDisplayActive = FALSE;
    return freedBytes;
}

static SbBool
mesh_int_field_contains(const SoMFInt32 &field, int value)
{
    for (int i = 0; i < field.getNum(); i++) {
	if (field[i] == value)
	    return TRUE;
    }
    return FALSE;
}

static float
mesh_distance_squared_to_segment(const SbVec3f &p, const SbVec3f &a, const SbVec3f &b)
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

    SbVec3f closest = a + ab * t;
    return (p - closest).sqrLength();
}

static int
mesh_nearest_face_vertex_slot(const SbVec3f &hit, const SbVec3f vertices[3])
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
mesh_nearest_face_edge_slot(const SbVec3f &hit, const SbVec3f vertices[3])
{
    static const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    int nearest = 0;
    float nearestDist = mesh_distance_squared_to_segment(hit,
			vertices[edges[0][0]], vertices[edges[0][1]]);
    for (int i = 1; i < 3; i++) {
	float dist = mesh_distance_squared_to_segment(hit,
		     vertices[edges[i][0]], vertices[edges[i][1]]);
	if (dist < nearestDist) {
	    nearest = i;
	    nearestDist = dist;
	}
    }
    return nearest;
}

SbBool
SoBRLMeshShape::isPrimitiveSelected(int primitiveIndex) const
{
    return this->selected.getValue() || mesh_int_field_contains(this->selectedPrimitive, primitiveIndex);
}

SbBool
SoBRLMeshShape::isPrimitiveHighlighted(int primitiveIndex) const
{
    return this->highlighted.getValue() || mesh_int_field_contains(this->highlightedPrimitive, primitiveIndex);
}

void
SoBRLMeshShape::setLodBackedMesh(SbBool enabled)
{
    this->lodBacked = enabled ? TRUE : FALSE;
}

SbBool
SoBRLMeshShape::isLodBackedMesh(void) const
{
    return this->lodBacked.getValue();
}

void
SoBRLMeshShape::setLodPreserveFullDetail(SbBool enabled)
{
    this->lodPreserveFullDetail = enabled ? TRUE : FALSE;
    if (!enabled)
	this->evictFullDetailMesh();
}

SbBool
SoBRLMeshShape::isLodPreserveFullDetailEnabled(void) const
{
    return this->lodPreserveFullDetail.getValue();
}

void
SoBRLMeshShape::setPickGeometryPolicy(PickGeometryPolicy policy)
{
    this->pickGeometryPolicy =
	(policy == PICK_FULL_DETAIL) ? PICK_FULL_DETAIL : PICK_DISPLAY_LEVEL;
}

SoBRLMeshShape::PickGeometryPolicy
SoBRLMeshShape::getPickGeometryPolicy(void) const
{
    return (this->pickGeometryPolicy.getValue() == PICK_FULL_DETAIL) ?
	   PICK_FULL_DETAIL : PICK_DISPLAY_LEVEL;
}

SbBool
SoBRLMeshShape::isLodDisplayActive(void) const
{
    return this->lodDisplayActive.getValue();
}

void
SoBRLMeshShape::makeLodRequest(BObolLodRequest &request,
			       const char *databaseId,
			       uint64_t databaseRevision,
			       uint64_t viewRevision,
			       uint64_t policyRevision,
			       int requestDrawMode,
			       const char *providerId,
			       const char *providerVersion,
			       int qualityTier) const
{
    request.clear();

    request.databaseId = databaseId ? databaseId : "";
    request.databaseRevision = databaseRevision;
    request.sourceRevision = this->sourceId.getValue();
    request.objectPath = this->sourcePath.getValue();
    request.objectName = this->sourceName.getValue();
    request.viewRevision = viewRevision;
    request.policyRevision = policyRevision;
    request.drawMode = requestDrawMode;
    request.lodPolicy = this->lodPolicy.getValue();
    request.providerId = providerId ? providerId : "";
    request.providerVersion = providerVersion ? providerVersion : "";
    request.qualityTier = qualityTier;
    const SoBRLMeshShape *geom = this->getSharedGeometrySource();
    const SbBool useFullDetail =
	(!geom->fullDetailPoint.empty() &&
	 !geom->fullDetailCoordIndex.empty()) ? TRUE : FALSE;
    const SbBool useSourceMetrics =
	(!useFullDetail && geom->sourceMeshMetricsValid) ? TRUE : FALSE;
    request.sourceCounts.faceCount = useFullDetail ?
				     static_cast<uint64_t>(geom->fullDetailCoordIndex.size() / 3) :
				     (useSourceMetrics ? geom->sourceMeshFaceCount :
				      static_cast<uint64_t>(geom->getTriangleCount()));
    request.sourceCounts.pointCount = useFullDetail ?
				      static_cast<uint64_t>(geom->fullDetailPoint.size()) :
				      (useSourceMetrics ? geom->sourceMeshPointCount :
				       static_cast<uint64_t>(geom->point.getNum()));

    SbBox3f bounds;
    bounds.makeEmpty();
    if (useFullDetail) {
	for (size_t i = 0; i < geom->fullDetailPoint.size(); i++)
	    bounds.extendBy(geom->fullDetailPoint[i]);
    } else if (useSourceMetrics) {
	bounds = geom->sourceMeshBounds;
    } else {
	for (int i = 0; i < geom->point.getNum(); i++)
	    bounds.extendBy(geom->point[i]);
    }
    request.bounds = bounds;
}

void
SoBRLMeshShape::clearStagedLodResult(void)
{
    SbBool keepDisplayPayload = FALSE;
    if (this->lodDisplayActive.getValue()) {
	if (this->hasFullDetailMesh())
	    this->restoreFullDetailMesh();
	else
	    keepDisplayPayload = TRUE;
    }
    this->lodDisplayActive = keepDisplayPayload;
    this->lodStagedAvailable = FALSE;
    this->lodResultKind = BOBOL_LOD_RESULT_NONE;
    this->lodQualityTier = BOBOL_LOD_QUALITY_METADATA;
    this->lodProviderStatus = BOBOL_LOD_PROVIDER_UNKNOWN;
    this->lodCacheKey = "";
    this->lodProviderId = "";
    this->lodProviderVersion = "";
    this->lodDiagnostic = "";
    this->lodDependencyPath.setNum(0);
    this->lodDependencyName.setNum(0);
    this->lodDependencySourceRevision.setNum(0);
    this->lodDependencySourceHash.setNum(0);
    this->lodDependencyQualityTier.setNum(0);
    this->lodDependencyOptional.setNum(0);
    this->lodAttributeName.setNum(0);
    this->lodAttributeValue.setNum(0);
    this->lodProxyKind = BOBOL_LOD_PROXY_NONE;
    this->lodProxyCenter = SbVec3f(0.0f, 0.0f, 0.0f);
    this->lodProxyAxisX = SbVec3f(1.0f, 0.0f, 0.0f);
    this->lodProxyAxisY = SbVec3f(0.0f, 1.0f, 0.0f);
    this->lodProxyAxisZ = SbVec3f(0.0f, 0.0f, 1.0f);
    this->lodProxyHalfExtents = SbVec3f(0.0f, 0.0f, 0.0f);
    this->lodProxyRadius = 0.0f;
}

SbBool
SoBRLMeshShape::applyStagedLodResult(const BObolLodResult &result,
				     const BObolLodRequest *expectedRequest)
{
    if (expectedRequest && !bobol_lod_result_matches_request(result,
	    *expectedRequest)) {
	this->clearStagedLodResult();
	this->lodProviderStatus = BOBOL_LOD_PROVIDER_STALE;
	this->lodDiagnostic = "stale staged LoD result rejected";
	return FALSE;
    }

    if (result.providerStatus == BOBOL_LOD_PROVIDER_STALE ||
	result.providerStatus == BOBOL_LOD_PROVIDER_CANCELLED ||
	result.providerStatus == BOBOL_LOD_PROVIDER_CACHE_MISS ||
	result.providerStatus == BOBOL_LOD_PROVIDER_ERROR) {
	this->clearStagedLodResult();
	this->lodResultKind = result.resultKind;
	this->lodQualityTier = result.qualityTier;
	this->lodProviderStatus = result.providerStatus;
	this->lodCacheKey = result.cacheKey.value;
	this->lodProviderId = result.request.providerId;
	this->lodProviderVersion = result.request.providerVersion;
	this->lodDiagnostic = result.diagnostic;
	return FALSE;
    }

    this->lodStagedAvailable = TRUE;
    this->lodResultKind = result.resultKind;
    this->lodQualityTier = result.qualityTier;
    this->lodProviderStatus = result.providerStatus;
    this->lodCacheKey = result.cacheKey.value;
    this->lodProviderId = result.request.providerId;
    this->lodProviderVersion = result.request.providerVersion;
    this->lodDiagnostic = result.diagnostic;

    this->lodFaceCount = mesh_lod_count_to_field(result.counts.faceCount);
    this->lodPointCount = mesh_lod_count_to_field(result.counts.pointCount);
    this->lodOriginalPointCount =
	mesh_lod_count_to_field(result.counts.originalPointCount);
    this->lodNormalCount = mesh_lod_count_to_field(result.counts.normalCount);
    this->lodHasSnappedPoints = result.hasSnappedPoints;
    this->lodHasNormals = result.hasNormals;
    if (!result.bounds.isEmpty()) {
	this->lodBoundsMin = result.bounds.getMin();
	this->lodBoundsMax = result.bounds.getMax();
    }
    if (result.resultKind == BOBOL_LOD_RESULT_MESH ||
	result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) {
	this->lodAvailable = TRUE;
	this->lodActiveCut = result.geometry.activeCut;
	if (result.mesh.isValid() &&
	    result.mesh.points.size() <=
	    static_cast<size_t>(std::numeric_limits<int>::max()) &&
	    result.mesh.coordIndex.size() <=
	    static_cast<size_t>(std::numeric_limits<int>::max())) {
	    if (result.resultKind == BOBOL_LOD_RESULT_MESH) {
		if (!this->sourceMeshMetricsValid)
		    this->updateSourceMeshMetricsFromFields();
		if (this->lodPreserveFullDetail.getValue())
		    this->captureFullDetailMesh();
		else
		    this->clearFullDetailMesh();
	    } else {
		this->clearFullDetailMesh();
	    }
	    const SbVec3f *normals = result.mesh.normals.empty() ?
				     NULL : result.mesh.normals.data();
	    const int normalCount = result.mesh.normals.empty() ? 0 :
				    static_cast<int>(result.mesh.normals.size());
	    this->setIndexedTriangleFields(
		result.mesh.points.data(),
		static_cast<int>(result.mesh.points.size()),
		result.mesh.coordIndex.data(),
		static_cast<int>(result.mesh.coordIndex.size()),
		normals,
		normalCount);
	    if (result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
		this->updateSourceMeshMetricsFromFields();
	    this->lodDisplayActive =
		(result.resultKind == BOBOL_LOD_RESULT_MESH) ? TRUE : FALSE;
	}
    }

    this->lodDependencyPath.setNum(0);
    this->lodDependencyName.setNum(0);
    this->lodDependencySourceRevision.setNum(0);
    this->lodDependencySourceHash.setNum(0);
    this->lodDependencyQualityTier.setNum(0);
    this->lodDependencyOptional.setNum(0);
    for (size_t i = 0; i < result.dependencies.size(); i++) {
	this->lodDependencyPath.set1Value(static_cast<int>(i),
					  result.dependencies[i].objectPath);
	this->lodDependencyName.set1Value(static_cast<int>(i),
					  result.dependencies[i].objectName);
	this->lodDependencySourceRevision.set1Value(static_cast<int>(i),
					 mesh_lod_uint64_string(result.dependencies[i].sourceRevision.value()));
	this->lodDependencySourceHash.set1Value(static_cast<int>(i),
						mesh_lod_uint64_string(result.dependencies[i].sourceContentHash));
	this->lodDependencyQualityTier.set1Value(static_cast<int>(i),
				      result.dependencies[i].requiredQualityTier);
	this->lodDependencyOptional.set1Value(static_cast<int>(i),
					      result.dependencies[i].optional);
    }

    this->lodAttributeName.setNum(0);
    this->lodAttributeValue.setNum(0);
    for (size_t i = 0; i < result.attributes.size(); i++) {
	this->lodAttributeName.set1Value(static_cast<int>(i),
					 result.attributes[i].name);
	this->lodAttributeValue.set1Value(static_cast<int>(i),
					  result.attributes[i].value);
    }

    this->lodProxyKind = result.proxy.kind;
    this->lodProxyCenter = result.proxy.center;
    this->lodProxyAxisX = result.proxy.axisX;
    this->lodProxyAxisY = result.proxy.axisY;
    this->lodProxyAxisZ = result.proxy.axisZ;
    this->lodProxyHalfExtents = result.proxy.halfExtents;
    this->lodProxyRadius = result.proxy.radius;

    return result.providerStatus == BOBOL_LOD_PROVIDER_ERROR ? FALSE : TRUE;
}

void
SoBRLMeshShape::updateSourceMeshMetricsFromFields(void)
{
    const SoBRLMeshShape *geom = this->getSharedGeometrySource();
    this->sourceMeshFaceCount = static_cast<uint64_t>(geom->getTriangleCount());
    this->sourceMeshPointCount = static_cast<uint64_t>(geom->point.getNum());
    this->sourceMeshBounds.makeEmpty();
    for (int i = 0; i < geom->point.getNum(); i++)
	this->sourceMeshBounds.extendBy(geom->point[i]);
    this->sourceMeshMetricsValid =
	(!this->sourceMeshBounds.isEmpty() &&
	 this->sourceMeshFaceCount > 0 &&
	 this->sourceMeshPointCount > 0) ? TRUE : FALSE;
}

void
SoBRLMeshShape::updateSourceMeshMetricsFromFullDetail(void)
{
    this->sourceMeshFaceCount =
	static_cast<uint64_t>(this->getFullDetailTriangleCount());
    this->sourceMeshPointCount =
	static_cast<uint64_t>(this->fullDetailPoint.size());
    this->sourceMeshBounds.makeEmpty();
    for (size_t i = 0; i < this->fullDetailPoint.size(); i++)
	this->sourceMeshBounds.extendBy(this->fullDetailPoint[i]);
    this->sourceMeshMetricsValid =
	(!this->sourceMeshBounds.isEmpty() &&
	 this->sourceMeshFaceCount > 0 &&
	 this->sourceMeshPointCount > 0) ? TRUE : FALSE;
}

void
SoBRLMeshShape::captureFullDetailMesh(void)
{
    if (this->hasFullDetailMesh() || this->lodDisplayActive.getValue())
	return;

    const SoBRLMeshShape *geom = this->getSharedGeometrySource();
    this->fullDetailPoint.clear();
    this->fullDetailCoordIndex.clear();
    if (geom->point.getNum() <= 0 || geom->coordIndex.getNum() <= 0)
	return;

    this->fullDetailPoint.reserve(static_cast<size_t>(geom->point.getNum()));
    for (int i = 0; i < geom->point.getNum(); i++)
	this->fullDetailPoint.push_back(geom->point[i]);

    this->fullDetailCoordIndex.reserve(
	static_cast<size_t>(geom->coordIndex.getNum()));
    for (int i = 0; i < geom->coordIndex.getNum(); i++)
	this->fullDetailCoordIndex.push_back(geom->coordIndex[i]);
}

void
SoBRLMeshShape::restoreFullDetailMesh(void)
{
    if (!this->hasFullDetailMesh())
	return;

    this->setIndexedTriangleFields(this->fullDetailPoint.data(),
				   static_cast<int>(this->fullDetailPoint.size()),
				   this->fullDetailCoordIndex.data(),
				   static_cast<int>(this->fullDetailCoordIndex.size()));
    this->updateSourceMeshMetricsFromFields();
    this->clearFullDetailMesh();
}

void
SoBRLMeshShape::clearFullDetailMesh(void)
{
    std::vector<SbVec3f>().swap(this->fullDetailPoint);
    std::vector<int32_t>().swap(this->fullDetailCoordIndex);
}

static SbBool
mesh_shape_gl_color(SoBRLMeshShape *shape, int primitiveIndex,
		    SbColor &color, float &alpha)
{
    SbBool overridden = TRUE;
    if (shape->isPrimitiveHighlighted(primitiveIndex)) {
	color = shape->highlightedColor.getValue();
	alpha = 1.0f;
    } else if (shape->isPrimitiveSelected(primitiveIndex)) {
	color = shape->selectedColor.getValue();
	alpha = 1.0f;
    } else if (shape->ghosted.getValue()) {
	color = shape->ghostedColor.getValue();
	alpha = 0.35f;
    } else if (shape->colorOverride.getValue()) {
	color = shape->color.getValue();
	alpha = 1.0f;
    } else {
	color = shape->materialColorValid.getValue() ?
	    shape->materialColor.getValue() : shape->color.getValue();
	alpha = 1.0f;
	overridden = FALSE;
    }

    const float transparency = std::max(0.0f,
	std::min(1.0f, shape->transparency.getValue()));
    alpha *= 1.0f - transparency;
    return overridden;
}

static void
mesh_shape_enable_transparency(SoBRLMeshShape *shape)
{
    const float transparency = std::max(0.0f,
	std::min(1.0f, shape->transparency.getValue()));
    if (transparency <= 0.0f && !shape->ghosted.getValue())
	return;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

static void
set_mesh_gl_material(SoBRLMeshShape *shape, int primitiveIndex)
{
    SbColor color;
    float alpha = 1.0f;
    (void)mesh_shape_gl_color(shape, primitiveIndex, color, alpha);

    const GLfloat ambient[4] = {
	color[0] * 0.2f,
	color[1] * 0.2f,
	color[2] * 0.2f,
	alpha
    };
    const GLfloat diffuse[4] = {
	color[0] * 0.6f,
	color[1] * 0.6f,
	color[2] * 0.6f,
	alpha
    };
    const GLfloat specular[4] = {
	color[0] * 0.2f,
	color[1] * 0.2f,
	color[2] * 0.2f,
	alpha
    };
    const GLfloat emission[4] = {0.0f, 0.0f, 0.0f, alpha};

    glColor4f(color[0], color[1], color[2], alpha);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
}

static void
set_mesh_gl_color(SoBRLMeshShape *shape, int primitiveIndex)
{
    SbColor color;
    float alpha = 1.0f;
    (void)mesh_shape_gl_color(shape, primitiveIndex, color, alpha);
    glColor4f(color[0], color[1], color[2], alpha);
}

static SbBool
mesh_shape_needs_per_primitive_color(const SoBRLMeshShape *shape)
{
    return shape &&
	   (shape->selectedPrimitive.getNum() > 0 ||
	    shape->highlightedPrimitive.getNum() > 0);
}

static void
mesh_shape_emit_triangles(SoBRLMeshShape *shape,
			  const BObolViewLodState::MeshPayload *viewPayload,
			  SbBool perPrimitiveColor)
{
    glBegin(GL_TRIANGLES);
    const int triangleCount = viewPayload ? viewPayload->getTriangleCount() :
			      shape->getTriangleCount();
    for (int i = 0; i < triangleCount; i++) {
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (viewPayload) {
	    if (!viewPayload->getTriangle(i, a, b, c))
		continue;
	} else if (!shape->getTriangle(i, a, b, c)) {
	    continue;
	}

	SbVec3f normalA;
	SbVec3f normalB;
	SbVec3f normalC;
	const int normalBase = i * 3;
	SbBool haveNormals = FALSE;
	if (viewPayload &&
		static_cast<size_t>(normalBase + 2) <
		viewPayload->mesh.normals.size()) {
	    normalA = viewPayload->mesh.normals[normalBase];
	    normalB = viewPayload->mesh.normals[normalBase + 1];
	    normalC = viewPayload->mesh.normals[normalBase + 2];
	    haveNormals = TRUE;
	} else if (!viewPayload) {
	    const SoBRLMeshShape *geom = shape->getGeometrySource();
	    if (normalBase + 2 < geom->normal.getNum()) {
		normalA = geom->normal[normalBase];
		normalB = geom->normal[normalBase + 1];
		normalC = geom->normal[normalBase + 2];
		haveNormals = TRUE;
	    }
	}

	if (!haveNormals) {
	    normalA = (b - a).cross(c - a);
	    if (normalA.length() > 0.0f)
		normalA.normalize();
	    else
		normalA = SbVec3f(0.0f, 0.0f, 1.0f);
	    normalB = normalA;
	    normalC = normalA;
	}

	if (perPrimitiveColor)
	    set_mesh_gl_color(shape, i);
	glNormal3f(normalA[0], normalA[1], normalA[2]);
	glVertex3f(a[0], a[1], a[2]);
	glNormal3f(normalB[0], normalB[1], normalB[2]);
	glVertex3f(b[0], b[1], b[2]);
	glNormal3f(normalC[0], normalC[1], normalC[2]);
	glVertex3f(c[0], c[1], c[2]);
    }
    glEnd();
}

static void
mesh_shape_emit_wire_lines(SoBRLMeshShape *shape,
			   const BObolViewLodState::MeshPayload *viewPayload,
			   SbBool perPrimitiveColor)
{
    glBegin(GL_LINES);
    const int triangleCount = viewPayload ? viewPayload->getTriangleCount() :
			      shape->getTriangleCount();
    for (int i = 0; i < triangleCount; i++) {
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (viewPayload) {
	    if (!viewPayload->getTriangle(i, a, b, c))
		continue;
	} else if (!shape->getTriangle(i, a, b, c)) {
	    continue;
	}

	if (perPrimitiveColor)
	    set_mesh_gl_color(shape, i);
	glVertex3f(a[0], a[1], a[2]);
	glVertex3f(b[0], b[1], b[2]);
	glVertex3f(b[0], b[1], b[2]);
	glVertex3f(c[0], c[1], c[2]);
	glVertex3f(c[0], c[1], c[2]);
	glVertex3f(a[0], a[1], a[2]);
    }
    glEnd();
}

/* Identify which geometry is currently baked into the cached display lists.
 * 'kind' distinguishes the triangle stream (shaded / hidden-line) from the
 * wire stream.  For a view-local LoD payload the cache key + active cut +
 * triangle count uniquely name the resident cut (the node fields do not
 * change, so notify() cannot see refinement).  For node-owned geometry the
 * triangle/normal counts suffice because notify() already discards the cache
 * whenever a field changes. */
static std::string
mesh_shape_geometry_signature(SoBRLMeshShape *shape, char kind,
			      const BObolViewLodState::MeshPayload *viewPayload)
{
    std::ostringstream out;
    out << kind << '|';
    if (viewPayload) {
	out << 'V' << '|'
	    << viewPayload->cacheKey.getString() << '|'
	    << viewPayload->activeCut << '|'
	    << viewPayload->getTriangleCount() << '|'
	    << static_cast<uint64_t>(viewPayload->mesh.normals.size());
    } else {
	const SoBRLMeshShape *geom = shape->getGeometrySource();
	out << 'N' << '|'
	    << shape->getTriangleCount() << '|'
	    << geom->normal.getNum();
    }
    return out.str();
}

/* Render the mesh geometry from a per-GL-context display list, building it on
 * first use (or after the baked signature changes) and re-calling it on every
 * subsequent frame.  This replaces per-frame immediate-mode submission of every
 * vertex, which dominated frame time for resident LoD meshes.  Only the raw
 * geometry (glBegin/glVertex/glNormal) is captured; material, lighting, line
 * width, polygon mode and transparency remain caller-set GL state applied when
 * the list is called, so the same triangle list serves both the shaded and
 * hidden-line passes. */
static void
mesh_shape_call_cached_geometry(SoBRLMeshShape *shape,
				SoGLRenderAction *action, char kind,
				const BObolViewLodState::MeshPayload *viewPayload)
{
    SoState *state = action->getState();
    const std::string signature =
	mesh_shape_geometry_signature(shape, kind, viewPayload);
    if (shape->renderListSignature != signature)
	shape->clearRenderLists();

    const int context = SoGLCacheContextElement::get(state);
    std::map<int, SoGLDisplayList *>::iterator found =
	shape->renderLists.find(context);
    if (found != shape->renderLists.end()) {
	found->second->call(state);
	return;
    }

    SoGLDisplayList *list =
	new SoGLDisplayList(state, SoGLDisplayList::DISPLAY_LIST);
    list->ref();
    shape->renderLists[context] = list;
    shape->renderListSignature = signature;
    list->open(state);
    if (kind == 'W')
	mesh_shape_emit_wire_lines(shape, viewPayload, FALSE);
    else
	mesh_shape_emit_triangles(shape, viewPayload, FALSE);
    list->close(state);
    list->call(state);
}

static void
mesh_shape_render_hidden_line(SoBRLMeshShape *shape,
			      SoGLRenderAction *action,
			      const BObolViewLodState::MeshPayload *viewPayload)
{
    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT |
		 GL_POLYGON_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT |
		 GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    mesh_shape_enable_transparency(shape);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthFunc(GL_LESS);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    const SbBool perPrimitiveColor =
	mesh_shape_needs_per_primitive_color(shape);
    if (perPrimitiveColor)
	mesh_shape_emit_triangles(shape, viewPayload, FALSE);
    else
	mesh_shape_call_cached_geometry(shape, action, 'T', viewPayload);

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDepthFunc(GL_LEQUAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    if (shape->lineWidth.getValue() > 0)
	glLineWidth(static_cast<GLfloat>(shape->lineWidth.getValue()));
    set_mesh_gl_material(shape, -1);
    if (perPrimitiveColor)
	mesh_shape_emit_triangles(shape, viewPayload, TRUE);
    else
	mesh_shape_call_cached_geometry(shape, action, 'T', viewPayload);
    glPopAttrib();
}

static void
mesh_shape_render_wire(SoBRLMeshShape *shape,
		       SoGLRenderAction *action,
		       const BObolViewLodState::MeshPayload *viewPayload)
{
    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT |
		 GL_LINE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    mesh_shape_enable_transparency(shape);
    if (shape->lineWidth.getValue() > 0)
	glLineWidth(static_cast<GLfloat>(shape->lineWidth.getValue()));
    set_mesh_gl_material(shape, -1);

    if (mesh_shape_needs_per_primitive_color(shape))
	mesh_shape_emit_wire_lines(shape, viewPayload, TRUE);
    else
	mesh_shape_call_cached_geometry(shape, action, 'W', viewPayload);
    glPopAttrib();
}

static SbBool
mesh_shape_proxy_corners(const BObolViewLodState::ProxyPayload *viewProxy,
			 SbVec3f corners[8])
{
    if (!viewProxy || !viewProxy->isValid())
	return FALSE;

    const BObolLodProxy &proxy = viewProxy->proxy;
    if (proxy.kind == BOBOL_LOD_PROXY_AABB) {
	if (proxy.bounds.isEmpty())
	    return FALSE;
	const SbVec3f bmin = proxy.bounds.getMin();
	const SbVec3f bmax = proxy.bounds.getMax();
	corners[0].setValue(bmin[0], bmin[1], bmin[2]);
	corners[1].setValue(bmax[0], bmin[1], bmin[2]);
	corners[2].setValue(bmax[0], bmax[1], bmin[2]);
	corners[3].setValue(bmin[0], bmax[1], bmin[2]);
	corners[4].setValue(bmin[0], bmin[1], bmax[2]);
	corners[5].setValue(bmax[0], bmin[1], bmax[2]);
	corners[6].setValue(bmax[0], bmax[1], bmax[2]);
	corners[7].setValue(bmin[0], bmax[1], bmax[2]);
	return TRUE;
    }

    if (proxy.kind == BOBOL_LOD_PROXY_OBB) {
	SbVec3f ax = proxy.axisX;
	SbVec3f ay = proxy.axisY;
	SbVec3f az = proxy.axisZ;
	if (ax.length() > 0.0f)
	    ax.normalize();
	if (ay.length() > 0.0f)
	    ay.normalize();
	if (az.length() > 0.0f)
	    az.normalize();
	ax *= proxy.halfExtents[0];
	ay *= proxy.halfExtents[1];
	az *= proxy.halfExtents[2];
	corners[0] = proxy.center - ax - ay - az;
	corners[1] = proxy.center + ax - ay - az;
	corners[2] = proxy.center + ax + ay - az;
	corners[3] = proxy.center - ax + ay - az;
	corners[4] = proxy.center - ax - ay + az;
	corners[5] = proxy.center + ax - ay + az;
	corners[6] = proxy.center + ax + ay + az;
	corners[7] = proxy.center - ax + ay + az;
	return TRUE;
    }

    return FALSE;
}

static void
mesh_shape_render_proxy(SoBRLMeshShape *shape,
			const BObolViewLodState::ProxyPayload *viewProxy)
{
    SbVec3f corners[8];
    if (!mesh_shape_proxy_corners(viewProxy, corners))
	return;

    static const int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT |
		 GL_LINE_BIT | GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    mesh_shape_enable_transparency(shape);
    if (shape->lineWidth.getValue() > 0)
	glLineWidth(static_cast<GLfloat>(shape->lineWidth.getValue()));
    set_mesh_gl_material(shape, -1);
    glBegin(GL_LINES);
    for (size_t i = 0; i < 12; i++) {
	const SbVec3f &a = corners[edges[i][0]];
	const SbVec3f &b = corners[edges[i][1]];
	glVertex3f(a[0], a[1], a[2]);
	glVertex3f(b[0], b[1], b[2]);
    }
    glEnd();
    glPopAttrib();
}

void
SoBRLMeshShape::GLRender(SoGLRenderAction *action)
{
    if (!this->visible.getValue() || !this->shouldGLRender(action))
	return;

    const BObolViewLodState::MeshPayload *viewPayload =
	bobol_view_lod_mesh_for_action(action, this);
    const BObolViewLodState::ProxyPayload *viewProxy =
	viewPayload ? NULL : bobol_view_lod_proxy_for_action(action, this);
    if (viewProxy) {
	mesh_shape_render_proxy(this, viewProxy);
	return;
    }

    if (this->hiddenLine.getValue() ||
	this->drawMode.getValue() == BOBOL_LOD_DRAW_HIDDEN_LINE) {
	mesh_shape_render_hidden_line(this, action, viewPayload);
	return;
    }

    if (this->drawMode.getValue() == BOBOL_LOD_DRAW_WIRE) {
	mesh_shape_render_wire(this, action, viewPayload);
	return;
    }

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT |
		 GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    mesh_shape_enable_transparency(this);
    set_mesh_gl_material(this, -1);
    const SbBool perPrimitiveColor =
	mesh_shape_needs_per_primitive_color(this);
    if (perPrimitiveColor) {
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	mesh_shape_emit_triangles(this, viewPayload, TRUE);
    } else {
	glDisable(GL_COLOR_MATERIAL);
	mesh_shape_call_cached_geometry(this, action, 'T', viewPayload);
    }
    glPopAttrib();
}

void
SoBRLMeshShape::computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center)
{
    box.makeEmpty();
    if (!this->visible.getValue()) {
	center = SbVec3f(0.0f, 0.0f, 0.0f);
	return;
    }

    const BObolViewLodState::MeshPayload *viewPayload =
	bobol_view_lod_mesh_for_action(action, this);
    if (viewPayload) {
	if (!viewPayload->bounds.isEmpty()) {
	    box = viewPayload->bounds;
	} else {
	    for (size_t i = 0; i < viewPayload->mesh.points.size(); i++)
		box.extendBy(viewPayload->mesh.points[i]);
	}
	center = box.isEmpty() ? SbVec3f(0.0f, 0.0f, 0.0f) : box.getCenter();
	return;
    }

    const BObolViewLodState::ProxyPayload *viewProxy =
	bobol_view_lod_proxy_for_action(action, this);
    if (viewProxy) {
	if (!viewProxy->bounds.isEmpty()) {
	    box = viewProxy->bounds;
	} else if (!viewProxy->proxy.bounds.isEmpty()) {
	    box = viewProxy->proxy.bounds;
	} else {
	    SbVec3f corners[8];
	    if (mesh_shape_proxy_corners(viewProxy, corners)) {
		for (size_t i = 0; i < 8; i++)
		    box.extendBy(corners[i]);
	    }
	}
	center = box.isEmpty() ? SbVec3f(0.0f, 0.0f, 0.0f) : box.getCenter();
	return;
    }

    const SoBRLMeshShape *geom = this->getGeometrySource();
    if (geom->point.getNum() > 0) {
	for (int i = 0; i < geom->point.getNum(); i++)
	    box.extendBy(geom->point[i]);
    } else if (geom->sourceMeshMetricsValid) {
	box = geom->sourceMeshBounds;
    }

    center = box.isEmpty() ? SbVec3f(0.0f, 0.0f, 0.0f) : box.getCenter();
}

void
SoBRLMeshShape::generatePrimitives(SoAction *action)
{
    if (!this->visible.getValue() || !this->selectable.getValue())
	return;

    const BObolViewLodState::MeshPayload *candidateViewPayload =
	bobol_view_lod_mesh_for_action(action, this);
    const BObolViewLodState::ProxyPayload *candidateViewProxy =
	candidateViewPayload ? NULL : bobol_view_lod_proxy_for_action(action,
	    this);
    const SbBool pickFullDetail =
	(this->getPickGeometryPolicy() == PICK_FULL_DETAIL &&
	 (this->lodDisplayActive.getValue() || candidateViewPayload ||
	  candidateViewProxy)) ?
	TRUE : FALSE;
    const SbBool usePreservedFullDetail =
	(pickFullDetail && this->hasFullDetailMesh()) ? TRUE : FALSE;
    if (pickFullDetail && !usePreservedFullDetail &&
	this->isLodBackedMesh())
	return;
    if (!pickFullDetail && candidateViewProxy)
	return;

    const BObolViewLodState::MeshPayload *viewPayload =
	pickFullDetail ? NULL : candidateViewPayload;
    const int triangleCount = usePreservedFullDetail ?
			      this->getFullDetailTriangleCount() :
			      (viewPayload ? viewPayload->getTriangleCount() :
			       this->getTriangleCount());
    const SoBRLMeshShape *geom = this->getGeometrySource();

    SoPrimitiveVertex v0;
    SoPrimitiveVertex v1;
    SoPrimitiveVertex v2;
    SoPointDetail p0;
    SoPointDetail p1;
    SoPointDetail p2;
    SoFaceDetail faceDetail;

    faceDetail.setNumPoints(3);

    v0.setDetail(&faceDetail);
    v1.setDetail(&faceDetail);
    v2.setDetail(&faceDetail);

    for (int i = 0; i < triangleCount; i++) {
	int ia = -1;
	int ib = -1;
	int ic = -1;
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (usePreservedFullDetail) {
	    if (!this->getFullDetailTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    a = geom->fullDetailPoint[static_cast<size_t>(ia)];
	    b = geom->fullDetailPoint[static_cast<size_t>(ib)];
	    c = geom->fullDetailPoint[static_cast<size_t>(ic)];
	} else if (viewPayload) {
	    if (!viewPayload->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    if (!viewPayload->getTriangle(i, a, b, c))
		continue;
	} else {
	    if (!this->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    a = geom->point[ia];
	    b = geom->point[ib];
	    c = geom->point[ic];
	}
	if (ia < 0 || ib < 0 || ic < 0)
	    continue;

	SbVec3f normalA;
	SbVec3f normalB;
	SbVec3f normalC;
	const int normalBase = i * 3;
	SbBool haveNormals = FALSE;
	if (viewPayload &&
		static_cast<size_t>(normalBase + 2) <
		viewPayload->mesh.normals.size()) {
	    normalA = viewPayload->mesh.normals[normalBase];
	    normalB = viewPayload->mesh.normals[normalBase + 1];
	    normalC = viewPayload->mesh.normals[normalBase + 2];
	    haveNormals = TRUE;
	} else if (!usePreservedFullDetail && !viewPayload) {
	    if (normalBase + 2 < geom->normal.getNum()) {
		normalA = geom->normal[normalBase];
		normalB = geom->normal[normalBase + 1];
		normalC = geom->normal[normalBase + 2];
		haveNormals = TRUE;
	    }
	}

	if (!haveNormals) {
	    normalA = (b - a).cross(c - a);
	    if (normalA.length() > 0.0f)
		normalA.normalize();
	    else
		normalA = SbVec3f(0.0f, 0.0f, 1.0f);
	    normalB = normalA;
	    normalC = normalA;
	}

	p0.setCoordinateIndex(ia);
	p1.setCoordinateIndex(ib);
	p2.setCoordinateIndex(ic);
	faceDetail.setPoint(0, &p0);
	faceDetail.setPoint(1, &p1);
	faceDetail.setPoint(2, &p2);
	faceDetail.setFaceIndex(i);
	v0.setPoint(a);
	v1.setPoint(b);
	v2.setPoint(c);
	v0.setNormal(normalA);
	v1.setNormal(normalB);
	v2.setNormal(normalC);
	this->invokeTriangleCallbacks(action, &v0, &v1, &v2);
    }
}

SoDetail *
SoBRLMeshShape::createTriangleDetail(SoRayPickAction *UNUSED(action),
				     const SoPrimitiveVertex *v1,
				     const SoPrimitiveVertex *v2,
				     const SoPrimitiveVertex *v3,
				     SoPickedPoint *pp)
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
    detail->setPrimitive(SoBRLPickDetail::FACE, -1);
    detail->setEditIntent(this->editIntentId.getValue(),
			  this->editIntentRole.getValue());

    const SoDetail *vertexDetail = v1 ? v1->getDetail() : NULL;
    if (vertexDetail && vertexDetail->isOfType(SoFaceDetail::getClassTypeId())) {
	const SoFaceDetail *faceDetail = static_cast<const SoFaceDetail *>(vertexDetail);
	detail->setPrimitive(SoBRLPickDetail::FACE, faceDetail->getFaceIndex());
	if (faceDetail->getNumPoints() >= 3) {
	    const SoPointDetail *p0 = faceDetail->getPoint(0);
	    const SoPointDetail *p1 = faceDetail->getPoint(1);
	    const SoPointDetail *p2 = faceDetail->getPoint(2);
	    int vertexIndices[3] = {
		p0 ? p0->getCoordinateIndex() : -1,
		p1 ? p1->getCoordinateIndex() : -1,
		p2 ? p2->getCoordinateIndex() : -1
	    };
	    detail->setFaceVertexIndices(vertexIndices[0], vertexIndices[1],
					 vertexIndices[2]);

	    if (v1 && v2 && v3) {
		static const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
		SbVec3f vertices[3] = {
		    v1->getPoint(),
		    v2->getPoint(),
		    v3->getPoint()
		};
		SbVec3f hit = pp ? pp->getObjectPoint(this) : v1->getPoint();
		int vertexSlot = mesh_nearest_face_vertex_slot(hit, vertices);
		int edgeSlot = mesh_nearest_face_edge_slot(hit, vertices);
		detail->setNearestFaceVertex(vertexSlot,
					     vertexIndices[vertexSlot]);
		detail->setNearestFaceEdge(edgeSlot,
					   vertexIndices[edges[edgeSlot][0]],
					   vertexIndices[edges[edgeSlot][1]]);
	    }
	}
    }
    if (pp)
	detail->setModelPoint(pp->getObjectPoint(this));
    else if (v1)
	detail->setModelPoint(v1->getPoint());
    return detail;
}
