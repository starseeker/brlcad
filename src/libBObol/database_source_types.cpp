/*        D A T A B A S E _ S O U R C E _ T Y P E S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

/** @file database_source_types.cpp
 *
 * Value-type defaults used by the database-source API.  Keeping these
 * mechanical constructors apart from source realization makes the ownership
 * and mutation implementation substantially easier to audit.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BViewQuery.h"

BObolDatabaseSourceSummary::BObolDatabaseSourceSummary(void) :
    valid(FALSE),
    path(""),
    instanceKey(""),
    parentInstanceKey(""),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION),
    displayName(""),
    hasParent(FALSE),
    drawTreeDepth(0),
    parentGroupPath(""),
    representationKey(""),
    representationMode(-1),
    drawMode(SoBRLDatabaseSource::WIREFRAME),
    sourceRevision(0),
    inputsRevision(0),
    viewRevision(0),
    realizedRevision(0),
    realizedSourceRevision(0),
    realizedInputsRevision(0),
    realizedViewRevision(0),
    realizationStatus(SoBRLDatabaseSource::UNREALIZED),
    realizationDiagnostic(""),
    realizationIdentity(""),
    realizationRoleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_NONE),
    realizationViewDependent(FALSE),
    realizationCsgLodEnabled(FALSE),
    realizationMeshLodEnabled(FALSE),
    realizationViewScale(0.0f),
    realizationLodScale(1.0f),
    realizationViewWidth(0),
    realizationViewHeight(0),
    realizationBotThreshold(0),
    realizationCurveScale(0.0f),
    realizationPointScale(0.0f),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    materialPolicy(SoBRLDatabaseSource::MATERIAL_INHERIT),
    databaseMetadataValid(FALSE),
    databaseRegionId(0),
    databaseAirCode(0),
    databaseMaterialId(0),
    databaseLos(0),
    databaseMaterialColorValid(FALSE),
    databaseMaterialColor(1.0f, 1.0f, 1.0f),
    databaseMaterialShader(""),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    selectedColor(1.0f, 1.0f, 1.0f),
    highlightedColor(1.0f, 1.0f, 0.0f),
    ghostedColor(0.55f, 0.55f, 0.55f),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    sourceBoundsValid(FALSE),
    sourceBoundsExact(FALSE),
    sourceBounds(),
    hasViewDependentCsgGeometry(FALSE),
    stale(TRUE),
    staleReason(SoBRLDatabaseSource::STALE_SOURCE),
    realizedShapeCount(0),
    realizedMeshCount(0),
    realizedMaterialObjectCount(0)
{
    this->sourceBounds.makeEmpty();
}

BObolAuxiliaryLineSetDisplayState::BObolAuxiliaryLineSetDisplayState(void) :
    valid(FALSE),
    drawMode(BOBOL_LOD_DRAW_WIRE),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0)
{
}

BObolCompactInstanceHandle::BObolCompactInstanceHandle(void) :
    sourceNodeId(0),
    instanceWord0(0),
    instanceWord1(0)
{
}

SbBool
BObolCompactInstanceHandle::isValid(void) const
{
    return sourceNodeId != 0 && (instanceWord0 != 0 || instanceWord1 != 0) ?
	TRUE : FALSE;
}

BObolCompactInstanceSummary::BObolCompactInstanceSummary(void) :
    valid(FALSE),
    sourceContentHash(0),
    sourceFaceCount(0),
    sourcePointCount(0),
    localToSource(SbMatrix::identity()),
    geometryIdentity(0),
    geometryRevision(0),
    appearanceRevision(0),
    placementRevision(0),
    visibilityRevision(0),
    selectionRevision(0),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialShader(""),
    appearanceColorValid(FALSE),
    appearanceColor(1.0f, 1.0f, 1.0f),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    wireGeometry(FALSE),
    pointGeometry(FALSE),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selectable(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    meshAssetBounds.makeEmpty();
}

BObolCompactLodInstanceSummary::BObolCompactLodInstanceSummary(void) :
    valid(FALSE),
    sourceContentHash(0),
    sourceFaceCount(0),
    sourcePointCount(0),
    brepSource(FALSE),
    meshAssetTessellationAbsTol(0.0),
    meshAssetTessellationRelTol(0.0),
    meshAssetTessellationNormTol(0.0),
    localToSource(SbMatrix::identity()),
    presentationLocalToSource(SbMatrix::identity()),
    presentationCornersValid(FALSE),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    meshAssetBounds.makeEmpty();
    localBounds.makeEmpty();
    presentationLocalBounds.makeEmpty();
}

BObolCompactLodProviderSummary::BObolCompactLodProviderSummary(void) :
    stagedSource(),
    lodAvailable(FALSE),
    lodActiveCut(-1),
    lodFaceCount(0),
    lodPointCount(0),
    lodHasNormals(FALSE)
{
}

BObolCompactLodPlanningSummary::BObolCompactLodPlanningSummary(void) :
    valid(FALSE),
    geometryRevision(0),
    placementRevision(0),
    sourceContentHash(0),
    sourceFaceCount(0),
    sourcePointCount(0),
    botSource(FALSE),
    brepSource(FALSE),
    meshAssetTessellationAbsTol(0.0),
    meshAssetTessellationRelTol(0.0),
    meshAssetTessellationNormTol(0.0),
    localToSource(SbMatrix::identity()),
    presentationLocalToSource(SbMatrix::identity()),
    presentationCornersValid(FALSE),
    meshGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE)
{
    localBounds.makeEmpty();
    presentationLocalBounds.makeEmpty();
}

BObolExternalLineSet::BObolExternalLineSet(void) :
    points(NULL),
    commands(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("line-set"),
    geometryKind("line")
{
}

BObolExternalPointSet::BObolExternalPointSet(void) :
    points(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("point-set"),
    geometryKind("point")
{
}

BObolExternalTriangleMesh::BObolExternalTriangleMesh(void) :
    points(NULL),
    pointCount(0),
    indices(NULL),
    indexCount(0),
    normals(NULL),
    normalCount(0),
    sourceType("indexed-face-set"),
    geometryKind("surface"),
    lodBacked(FALSE)
{
}

BObolExternalAnnotationSegment::BObolExternalAnnotationSegment(void) :
    kind(SEGMENT_NONE),
    lineStart(0),
    lineEnd(0),
    textRefPoint(0),
    text(NULL)
{
}

BObolExternalAnnotation::BObolExternalAnnotation(void) :
    basePoint(0.0f, 0.0f, 0.0f),
    linePoints(NULL),
    lineCommands(NULL),
    preciseLinePoints(NULL),
    linePointCount(0),
    annotationPoints(NULL),
    preciseAnnotationPoints(NULL),
    annotationPointCount(0),
    segments(NULL),
    segmentCount(0),
    sourceType("annotation"),
    geometryKind("annotation")
{
}

BObolDatabaseSourceDisplayPatch::BObolDatabaseSourceDisplayPatch(void) :
    visibleValid(FALSE),
    visible(TRUE),
    selectedValid(FALSE),
    selected(FALSE),
    highlightedValid(FALSE),
    highlighted(FALSE),
    lineStyleValid(FALSE),
    lineStyle(0),
    lineWidthValid(FALSE),
    lineWidth(0),
    transparencyValid(FALSE),
    transparency(0.0f),
    colorOverrideValid(FALSE),
    colorOverride(FALSE),
    colorValid(FALSE),
    color(1.0f, 1.0f, 1.0f),
    selectedColorValid(FALSE),
    selectedColor(1.0f, 1.0f, 1.0f),
    highlightedColorValid(FALSE),
    highlightedColor(1.0f, 1.0f, 0.0f),
    ghostedColorValid(FALSE),
    ghostedColor(0.55f, 0.55f, 0.55f)
{
}

BObolRealizedShapeSummary::BObolRealizedShapeSummary(void) :
    valid(FALSE),
    shapeKind(SHAPE_UNKNOWN),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    displayName(""),
    geometryName(""),
    cacheIdentity(""),
    sourceIdentity(""),
    databaseIntent(FALSE),
    overlayIntent(FALSE),
    hudIntent(FALSE),
    localSource(FALSE),
    sharedSource(FALSE),
    nonDatabaseSource(FALSE),
    drawMode(0),
    recordRole(""),
    geometryKind(""),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialShader(""),
    materialRevision(0),
    visible(FALSE),
    selectable(FALSE),
    selected(FALSE),
    highlighted(FALSE),
    ghosted(FALSE),
    hiddenLine(FALSE),
    editEmphasis(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    editIntentId(""),
    editIntentRole(""),
    lodPolicy(0),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    lodAvailable(FALSE),
    lodActiveCut(-1),
    lodFaceCount(0),
    lodPointCount(0),
    lodOriginalPointCount(0),
    lodNormalCount(0),
    lodHasSnappedPoints(FALSE),
    lodHasNormals(FALSE),
    lodBoundsMin(0.0f, 0.0f, 0.0f),
    lodBoundsMax(0.0f, 0.0f, 0.0f),
    pointCount(0),
    commandCount(0),
    segmentCount(0),
    pointPrimitiveCount(0),
    triangleCount(0),
    indexCount(0),
    boundsValid(FALSE),
    bounds()
{
    this->bounds.makeEmpty();
}

BObolCompactOccurrence::BObolCompactOccurrence(void) :
    geometry(),
    summary(),
    geometryTransform(SbMatrix::identity()),
    localTransform(SbMatrix::identity()),
    viewDependentCsgGeometry(FALSE),
    lodBacked(FALSE),
    sourceMeshRequestValid(FALSE),
    sourceMeshRequest(),
    occurrenceIndex(0),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION)
{
}

BObolCompactSourceProfile::BObolCompactSourceProfile(void) :
    valid(FALSE), occurrenceCount(0), uniqueAssetCount(0),
    encodedSourceBytes(0), largestAssetBytes(0), reusedOccurrenceCount(0)
{
}

BObolRealizedMaterialSummary::BObolRealizedMaterialSummary(void) :
    valid(FALSE),
    sourcePath(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    materialName(""),
    parentName(""),
    materialSource(""),
    propertyCount(0)
{
}

BObolSceneTreeSummary::BObolSceneTreeSummary(void) :
    valid(FALSE),
    nodeKind(NODE_UNKNOWN),
    isGroup(FALSE),
    isShape(FALSE),
    isDatabaseSource(FALSE),
    isMaterialObject(FALSE),
    hasParent(FALSE),
    drawTreeDepth(0),
    childCount(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    displayName(""),
    geometryName("")
{
}

BObolSceneDisplaySummary::BObolSceneDisplaySummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    isDatabaseSource(FALSE),
    hasDrawIntent(FALSE),
    intentPath(""),
    intentDrawMode(-1),
    visible(TRUE),
    selected(FALSE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0),
    drawMode(0),
    materialValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BObolSceneMaterialSummary::BObolSceneMaterialSummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    materialValid(FALSE),
    materialRevision(0),
    materialColor(1.0f, 1.0f, 1.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BObolSceneBoundsSummary::BObolSceneBoundsSummary(void) :
    valid(FALSE),
    nodeKind(BObolSceneTreeSummary::NODE_UNKNOWN),
    boundsValid(FALSE),
    bounds(),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
    this->bounds.makeEmpty();
}
