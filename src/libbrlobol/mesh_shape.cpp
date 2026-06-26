/*                   M E S H _ S H A P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/lod_realization.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/pick_detail.h"

#include <Inventor/SbVec3f.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoFaceDetail.h>
#include <Inventor/details/SoPointDetail.h>
#include <Inventor/gl.h>

#include <limits>
#include <sstream>
#include <stdint.h>

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
    SO_NODE_ADD_FIELD(sourcePath, (""));
    SO_NODE_ADD_FIELD(sourceName, (""));
    SO_NODE_ADD_FIELD(sourceType, (""));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(displayName, (""));
    SO_NODE_ADD_FIELD(geometryName, (""));
    SO_NODE_ADD_FIELD(cacheIdentity, (""));
    SO_NODE_ADD_FIELD(sourceIdentity, (""));
    SO_NODE_ADD_FIELD(ownerSourcePath, (""));
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
    SO_NODE_ADD_FIELD(selectedColor, (SbColor(0.0f, 0.75f, 1.0f)));
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
    SO_NODE_ADD_FIELD(lodActiveLevel, (-1));
    SO_NODE_ADD_FIELD(lodFaceCount, (0));
    SO_NODE_ADD_FIELD(lodPointCount, (0));
    SO_NODE_ADD_FIELD(lodOriginalPointCount, (0));
    SO_NODE_ADD_FIELD(lodNormalCount, (0));
    SO_NODE_ADD_FIELD(lodHasSnappedPoints, (FALSE));
    SO_NODE_ADD_FIELD(lodHasNormals, (FALSE));
    SO_NODE_ADD_FIELD(lodBoundsMin, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodBoundsMax, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(lodStagedAvailable, (FALSE));
    SO_NODE_ADD_FIELD(lodResultKind, (BRLOBOL_LOD_RESULT_NONE));
    SO_NODE_ADD_FIELD(lodQualityTier, (BRLOBOL_LOD_QUALITY_METADATA));
    SO_NODE_ADD_FIELD(lodProviderStatus, (BRLOBOL_LOD_PROVIDER_UNKNOWN));
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
    SO_NODE_ADD_FIELD(lodProxyKind, (BRLOBOL_LOD_PROXY_NONE));
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
}

SoBRLMeshShape::~SoBRLMeshShape(void)
{
}

void
SoBRLMeshShape::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLMeshShape, SoShape, "Shape");
}

void
SoBRLMeshShape::setIndexedTriangles(const SbVec3f *points, int pointCount,
	const int32_t *indices, int indexCount)
{
    this->clearFullDetailMesh();
    this->lodDisplayActive = FALSE;
    this->setIndexedTriangleFields(points, pointCount, indices, indexCount);
    this->updateSourceMeshMetricsFromFields();
}

void
SoBRLMeshShape::setIndexedTriangleFields(const SbVec3f *points, int pointCount,
	const int32_t *indices, int indexCount)
{
    this->point.setNum(0);
    this->coordIndex.setNum(0);
    if (!points || pointCount <= 0 || !indices || indexCount <= 0)
	return;

    this->point.setValues(0, pointCount, points);
    this->coordIndex.setValues(0, indexCount, indices);
}

int
SoBRLMeshShape::getTriangleCount(void) const
{
    return this->coordIndex.getNum() / 3;
}

SbBool
SoBRLMeshShape::getTriangle(int triangleIndex, SbVec3f &a, SbVec3f &b, SbVec3f &c) const
{
    int ia = -1;
    int ib = -1;
    int ic = -1;
    if (!this->getTriangleVertexIndices(triangleIndex, ia, ib, ic))
	return FALSE;

    a = this->point[ia];
    b = this->point[ib];
    c = this->point[ic];
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

    int base = triangleIndex * 3;
    if (base + 2 >= this->coordIndex.getNum())
	return FALSE;

    indexA = this->coordIndex[base];
    indexB = this->coordIndex[base + 1];
    indexC = this->coordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	    indexA >= this->point.getNum() ||
	    indexB >= this->point.getNum() ||
	    indexC >= this->point.getNum()) {
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
    return (!this->fullDetailPoint.empty() &&
	    !this->fullDetailCoordIndex.empty()) ? TRUE : FALSE;
}

int
SoBRLMeshShape::getFullDetailTriangleCount(void) const
{
    return static_cast<int>(this->fullDetailCoordIndex.size() / 3);
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

    a = this->fullDetailPoint[static_cast<size_t>(ia)];
    b = this->fullDetailPoint[static_cast<size_t>(ib)];
    c = this->fullDetailPoint[static_cast<size_t>(ic)];
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

    size_t base = static_cast<size_t>(triangleIndex) * 3;
    if (base + 2 >= this->fullDetailCoordIndex.size())
	return FALSE;

    indexA = this->fullDetailCoordIndex[base];
    indexB = this->fullDetailCoordIndex[base + 1];
    indexC = this->fullDetailCoordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	    static_cast<size_t>(indexA) >= this->fullDetailPoint.size() ||
	    static_cast<size_t>(indexB) >= this->fullDetailPoint.size() ||
	    static_cast<size_t>(indexC) >= this->fullDetailPoint.size()) {
	indexA = -1;
	indexB = -1;
	indexC = -1;
	return FALSE;
    }

    return TRUE;
}

SbBool
SoBRLMeshShape::makeSourceMeshRequest(BRLObolSourceMeshRequest &request) const
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
    request.lodActiveLevel = this->lodActiveLevel.getValue();
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

    if (this->hasFullDetailMesh()) {
	request.faceCount = static_cast<uint64_t>(this->getFullDetailTriangleCount());
	request.pointCount = static_cast<uint64_t>(this->fullDetailPoint.size());
	for (size_t i = 0; i < this->fullDetailPoint.size(); i++)
	    request.bounds.extendBy(this->fullDetailPoint[i]);
    } else if (this->sourceMeshMetricsValid) {
	request.faceCount = this->sourceMeshFaceCount;
	request.pointCount = this->sourceMeshPointCount;
	request.bounds = this->sourceMeshBounds;
    } else {
	request.faceCount = static_cast<uint64_t>(this->getTriangleCount());
	request.pointCount = static_cast<uint64_t>(this->point.getNum());
	for (int i = 0; i < this->point.getNum(); i++)
	    request.bounds.extendBy(this->point[i]);
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
		this->lodResultKind.getValue() == BRLOBOL_LOD_RESULT_MESH))
	return FALSE;

    BRLObolSourceMeshRequest request;
    return this->makeSourceMeshRequest(request);
}

size_t
SoBRLMeshShape::estimateDisplayMeshBytes(void) const
{
    size_t bytes = static_cast<size_t>(this->point.getNum()) * sizeof(SbVec3f);
    bytes += static_cast<size_t>(this->coordIndex.getNum()) * sizeof(int32_t);
    return bytes;
}

size_t
SoBRLMeshShape::estimateFullDetailMeshBytes(void) const
{
    size_t bytes = this->fullDetailPoint.capacity() * sizeof(SbVec3f);
    bytes += this->fullDetailCoordIndex.capacity() * sizeof(int32_t);
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
    size_t freedBytes = this->estimateFullDetailMeshBytes();
    if (freedBytes == 0)
	return 0;

    if (this->hasFullDetailMesh())
	this->updateSourceMeshMetricsFromFullDetail();
    this->clearFullDetailMesh();
    return freedBytes;
}

size_t
SoBRLMeshShape::evictActiveDisplayMesh(void)
{
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
SoBRLMeshShape::makeLodRequest(BRLObolLodRequest &request,
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
    const SbBool useFullDetail = this->hasFullDetailMesh();
    const SbBool useSourceMetrics =
	(!useFullDetail && this->sourceMeshMetricsValid) ? TRUE : FALSE;
    request.sourceCounts.faceCount = useFullDetail ?
	static_cast<uint64_t>(this->getFullDetailTriangleCount()) :
	(useSourceMetrics ? this->sourceMeshFaceCount :
	    static_cast<uint64_t>(this->getTriangleCount()));
    request.sourceCounts.pointCount = useFullDetail ?
	static_cast<uint64_t>(this->fullDetailPoint.size()) :
	(useSourceMetrics ? this->sourceMeshPointCount :
	    static_cast<uint64_t>(this->point.getNum()));

    SbBox3f bounds;
    bounds.makeEmpty();
    if (useFullDetail) {
	for (size_t i = 0; i < this->fullDetailPoint.size(); i++)
	    bounds.extendBy(this->fullDetailPoint[i]);
    } else if (useSourceMetrics) {
	bounds = this->sourceMeshBounds;
    } else {
	for (int i = 0; i < this->point.getNum(); i++)
	    bounds.extendBy(this->point[i]);
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
    this->lodResultKind = BRLOBOL_LOD_RESULT_NONE;
    this->lodQualityTier = BRLOBOL_LOD_QUALITY_METADATA;
    this->lodProviderStatus = BRLOBOL_LOD_PROVIDER_UNKNOWN;
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
    this->lodProxyKind = BRLOBOL_LOD_PROXY_NONE;
    this->lodProxyCenter = SbVec3f(0.0f, 0.0f, 0.0f);
    this->lodProxyAxisX = SbVec3f(1.0f, 0.0f, 0.0f);
    this->lodProxyAxisY = SbVec3f(0.0f, 1.0f, 0.0f);
    this->lodProxyAxisZ = SbVec3f(0.0f, 0.0f, 1.0f);
    this->lodProxyHalfExtents = SbVec3f(0.0f, 0.0f, 0.0f);
    this->lodProxyRadius = 0.0f;
}

SbBool
SoBRLMeshShape::applyStagedLodResult(const BRLObolLodResult &result,
	const BRLObolLodRequest *expectedRequest)
{
    if (expectedRequest && !brlobol_lod_result_matches_request(result,
	    *expectedRequest)) {
	this->clearStagedLodResult();
	this->lodProviderStatus = BRLOBOL_LOD_PROVIDER_STALE;
	this->lodDiagnostic = "stale staged LoD result rejected";
	return FALSE;
    }

    if (result.providerStatus == BRLOBOL_LOD_PROVIDER_STALE ||
	    result.providerStatus == BRLOBOL_LOD_PROVIDER_CANCELLED ||
	    result.providerStatus == BRLOBOL_LOD_PROVIDER_CACHE_MISS ||
	    result.providerStatus == BRLOBOL_LOD_PROVIDER_ERROR) {
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
    if (result.resultKind == BRLOBOL_LOD_RESULT_MESH ||
	    result.resultKind == BRLOBOL_LOD_RESULT_FULL_DETAIL) {
	this->lodAvailable = TRUE;
	this->lodActiveLevel = result.geometry.activeLevel;
	if (result.mesh.isValid() &&
		result.mesh.points.size() <=
		static_cast<size_t>(std::numeric_limits<int>::max()) &&
		result.mesh.coordIndex.size() <=
		static_cast<size_t>(std::numeric_limits<int>::max())) {
	    if (result.resultKind == BRLOBOL_LOD_RESULT_MESH) {
		if (!this->sourceMeshMetricsValid)
		    this->updateSourceMeshMetricsFromFields();
		if (this->lodPreserveFullDetail.getValue())
		    this->captureFullDetailMesh();
		else
		    this->clearFullDetailMesh();
	    } else {
		this->clearFullDetailMesh();
	    }
	    this->setIndexedTriangleFields(result.mesh.points.data(),
		    static_cast<int>(result.mesh.points.size()),
		    result.mesh.coordIndex.data(),
		    static_cast<int>(result.mesh.coordIndex.size()));
	    if (result.resultKind == BRLOBOL_LOD_RESULT_FULL_DETAIL)
		this->updateSourceMeshMetricsFromFields();
	    this->lodDisplayActive =
		(result.resultKind == BRLOBOL_LOD_RESULT_MESH) ? TRUE : FALSE;
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
		mesh_lod_uint64_string(result.dependencies[i].sourceRevision));
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

    return result.providerStatus == BRLOBOL_LOD_PROVIDER_ERROR ? FALSE : TRUE;
}

void
SoBRLMeshShape::updateSourceMeshMetricsFromFields(void)
{
    this->sourceMeshFaceCount = static_cast<uint64_t>(this->getTriangleCount());
    this->sourceMeshPointCount = static_cast<uint64_t>(this->point.getNum());
    this->sourceMeshBounds.makeEmpty();
    for (int i = 0; i < this->point.getNum(); i++)
	this->sourceMeshBounds.extendBy(this->point[i]);
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

    this->fullDetailPoint.clear();
    this->fullDetailCoordIndex.clear();
    if (this->point.getNum() <= 0 || this->coordIndex.getNum() <= 0)
	return;

    this->fullDetailPoint.reserve(static_cast<size_t>(this->point.getNum()));
    for (int i = 0; i < this->point.getNum(); i++)
	this->fullDetailPoint.push_back(this->point[i]);

    this->fullDetailCoordIndex.reserve(
	    static_cast<size_t>(this->coordIndex.getNum()));
    for (int i = 0; i < this->coordIndex.getNum(); i++)
	this->fullDetailCoordIndex.push_back(this->coordIndex[i]);
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

static void
set_mesh_gl_color(SoBRLMeshShape *shape, int primitiveIndex)
{
    if (shape->isPrimitiveHighlighted(primitiveIndex)) {
	const SbColor &c = shape->highlightedColor.getValue();
	glColor3f(c[0], c[1], c[2]);
    } else if (shape->isPrimitiveSelected(primitiveIndex)) {
	const SbColor &c = shape->selectedColor.getValue();
	glColor3f(c[0], c[1], c[2]);
    } else if (shape->ghosted.getValue()) {
	const SbColor &c = shape->ghostedColor.getValue();
	glColor4f(c[0], c[1], c[2], 0.35f);
    } else if (shape->colorOverride.getValue()) {
	const SbColor &c = shape->color.getValue();
	glColor3f(c[0], c[1], c[2]);
    }
}

void
SoBRLMeshShape::GLRender(SoGLRenderAction *action)
{
    if (!this->visible.getValue() || !this->shouldGLRender(action))
	return;

    glPushAttrib(GL_CURRENT_BIT);
    set_mesh_gl_color(this, -1);

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < this->getTriangleCount(); i++) {
	SbVec3f a;
	SbVec3f b;
	SbVec3f c;
	if (!this->getTriangle(i, a, b, c))
	    continue;

	SbVec3f normal = (b - a).cross(c - a);
	if (normal.length() > 0.0f)
	    normal.normalize();
	else
	    normal = SbVec3f(0.0f, 0.0f, 1.0f);

	set_mesh_gl_color(this, i);
	glNormal3f(normal[0], normal[1], normal[2]);
	glVertex3f(a[0], a[1], a[2]);
	glVertex3f(b[0], b[1], b[2]);
	glVertex3f(c[0], c[1], c[2]);
    }
    glEnd();
    glPopAttrib();
}

void
SoBRLMeshShape::computeBBox(SoAction *UNUSED(action), SbBox3f &box, SbVec3f &center)
{
    box.makeEmpty();
    if (!this->visible.getValue()) {
	center = SbVec3f(0.0f, 0.0f, 0.0f);
	return;
    }

    if (this->point.getNum() > 0) {
	for (int i = 0; i < this->point.getNum(); i++)
	    box.extendBy(this->point[i]);
    } else if (this->sourceMeshMetricsValid) {
	box = this->sourceMeshBounds;
    }

    center = box.isEmpty() ? SbVec3f(0.0f, 0.0f, 0.0f) : box.getCenter();
}

void
SoBRLMeshShape::generatePrimitives(SoAction *action)
{
    if (!this->visible.getValue() || !this->selectable.getValue())
	return;

    const SbBool pickFullDetail =
	(this->getPickGeometryPolicy() == PICK_FULL_DETAIL &&
	 this->lodDisplayActive.getValue()) ? TRUE : FALSE;
    if (pickFullDetail && !this->hasFullDetailMesh())
	return;

    const int triangleCount = pickFullDetail ?
	this->getFullDetailTriangleCount() : this->getTriangleCount();

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
	if (pickFullDetail) {
	    if (!this->getFullDetailTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    a = this->fullDetailPoint[static_cast<size_t>(ia)];
	    b = this->fullDetailPoint[static_cast<size_t>(ib)];
	    c = this->fullDetailPoint[static_cast<size_t>(ic)];
	} else {
	    if (!this->getTriangleVertexIndices(i, ia, ib, ic))
		continue;
	    a = this->point[ia];
	    b = this->point[ib];
	    c = this->point[ic];
	}
	if (ia < 0 || ib < 0 || ic < 0)
	    continue;

	SbVec3f normal = (b - a).cross(c - a);
	if (normal.length() > 0.0f)
	    normal.normalize();
	else
	    normal = SbVec3f(0.0f, 0.0f, 1.0f);

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
	v0.setNormal(normal);
	v1.setNormal(normal);
	v2.setNormal(normal);
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
