/*             B D A T A B A S E S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BDatabaseSource.h */

#ifndef BOBOL_BDATABASESOURCE_H
#define BOBOL_BDATABASESOURCE_H

#include "BObol/BDefines.h"
#include "BObol/BSourceMeshRequest.h"
#include "BObol/BViewLod.h"

#include <stdint.h>

#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

#include <memory>
#include <vector>

class SoBRLVListShape;
class SoBRLMeshShape;
class SoBRLMaterialObject;
class SoBRLCadAssembly;
class SoBRLCadRenderBatch;
class SoBRLExportAction;
class SoBRLMeasureAction;
class SoBRLSnapAction;
class SoFieldSensor;
class SoCallbackAction;
class SoGetBoundingBoxAction;
class SoGLRenderAction;
class SoRayPickAction;
class SoSensor;
struct db_i;
struct BObolMeshLod;
struct db_full_path;
struct rt_db_internal;
struct bg_tess_tol;
struct bn_tol;
struct BObolDatabaseSourceRealizationCache;
struct BObolCompactInstanceIndex;
struct BObolCompactInstanceEntry;
struct BObolCadBatchBuildState;
namespace Obol { struct PartGeometry; }

BOBOL_EXPORT SbBool bobol_database_source_fullpath_material_color(
	struct db_i *dbip,
	const struct db_full_path *pathp,
	SbColor &color);
BOBOL_EXPORT SbBool bobol_database_source_path_material_color(
	struct db_i *dbip,
	const char *path,
	SbColor &color);

struct BOBOL_EXPORT BObolDatabaseSourceSummary {
    BObolDatabaseSourceSummary(void);

    SbBool valid;
    SbString path;
    SbString instanceKey;
    SbString parentInstanceKey;
    uint32_t occurrenceIndex;
    int booleanOperation;
    SbString displayName;
    SbBool hasParent;
    int drawTreeDepth;
    SbString parentGroupPath;
    SbString representationKey;
    int representationMode;
    int drawMode;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    uint32_t viewRevision;
    uint32_t realizedRevision;
    uint32_t realizedSourceRevision;
    uint32_t realizedInputsRevision;
    uint32_t realizedViewRevision;
    int realizationStatus;
    SbString realizationDiagnostic;
    SbString realizationIdentity;
    int realizationRoleFlags;
    SbBool realizationViewDependent;
    SbBool realizationCsgLodEnabled;
    SbBool realizationMeshLodEnabled;
    float realizationViewScale;
    float realizationLodScale;
    int realizationViewWidth;
    int realizationViewHeight;
    uint32_t realizationBotThreshold;
    float realizationCurveScale;
    float realizationPointScale;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
    int materialPolicy;
    SbBool databaseMetadataValid;
    int databaseRegionId;
    int databaseAirCode;
    int databaseMaterialId;
    int databaseLos;
    SbBool databaseMaterialColorValid;
    SbColor databaseMaterialColor;
    SbString databaseMaterialShader;
    SbBool colorOverride;
    SbColor color;
    SbColor selectedColor;
    SbColor highlightedColor;
    SbColor ghostedColor;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
    SbBool sourceBoundsValid;
    SbBox3f sourceBounds;
    SbBool stale;
    uint32_t staleReason;
    int realizedShapeCount;
    int realizedMeshCount;
    int realizedMaterialObjectCount;
};

struct BOBOL_EXPORT BObolCompactInstanceHandle {
    BObolCompactInstanceHandle(void);

    uint64_t sourceNodeId;
    uint64_t instanceWord0;
    uint64_t instanceWord1;

    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolCompactInstanceSummary {
    BObolCompactInstanceSummary(void);

    SbBool valid;
    BObolCompactInstanceHandle handle;
    SbString path;
    SbString sourceName;
    SbString sourceInstanceKey;
    SbMatrix localToSource;
    SbBox3f localBounds;
    uint64_t geometryIdentity;
    uint64_t geometryRevision;
    uint64_t appearanceRevision;
    uint64_t placementRevision;
    uint64_t visibilityRevision;
    uint64_t selectionRevision;
    uint32_t occurrenceIndex;
    int booleanOperation;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    SbColor materialColor;
    SbString materialShader;
    /* Effective compact-instance presentation after selected/highlighted
     * state has been applied. */
    SbBool appearanceColorValid;
    SbColor appearanceColor;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool wireGeometry;
    SbBool pointGeometry;
    SbBool meshGeometry;
    SbBool lodBacked;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool highlighted;
};

struct BOBOL_EXPORT BObolAuxiliaryLineSetDisplayState {
    BObolAuxiliaryLineSetDisplayState(void);

    SbBool valid;
    int drawMode;
    SbBool visible;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
};

/*
 * External primary geometry published to a database source is source-local.
 * Callers must not pre-apply the source drawMatrix to these point arrays.
 * Bridges from legacy instance/world-space producers are responsible for
 * normalizing points back to source-local coordinates before publication.
 */
struct BOBOL_EXPORT BObolExternalLineSet {
    BObolExternalLineSet(void);

    const SbVec3f *points;
    const int32_t *commands;
    const double *precisePoints;
    int count;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolExternalPointSet {
    BObolExternalPointSet(void);

    const SbVec3f *points;
    const double *precisePoints;
    int count;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolExternalTriangleMesh {
    BObolExternalTriangleMesh(void);

    const SbVec3f *points;
    int pointCount;
    const int32_t *indices;
    int indexCount;
    /* Optional triangle-corner normals.  When present normalCount must match
     * indexCount (three normals per triangle). */
    const SbVec3f *normals;
    int normalCount;
    const char *sourceType;
    const char *geometryKind;
    SbBool lodBacked;
};

struct BOBOL_EXPORT BObolExternalAnnotationSegment {
    enum SegmentKind {
	SEGMENT_NONE = 0,
	SEGMENT_LINE = 1,
	SEGMENT_TEXT = 2
    };

    BObolExternalAnnotationSegment(void);

    int kind;
    int lineStart;
    int lineEnd;
    int textRefPoint;
    const char *text;
};

struct BOBOL_EXPORT BObolExternalAnnotation {
    BObolExternalAnnotation(void);

    SbVec3f basePoint;
    const SbVec3f *linePoints;
    const int32_t *lineCommands;
    const double *preciseLinePoints;
    int linePointCount;
    const SbVec3f *annotationPoints;
    const double *preciseAnnotationPoints;
    int annotationPointCount;
    const BObolExternalAnnotationSegment *segments;
    int segmentCount;
    const char *sourceType;
    const char *geometryKind;
};

struct BOBOL_EXPORT BObolDatabaseSourceDisplayPatch {
    BObolDatabaseSourceDisplayPatch(void);

    SbBool visibleValid;
    SbBool visible;
    SbBool selectedValid;
    SbBool selected;
    SbBool highlightedValid;
    SbBool highlighted;
    SbBool lineStyleValid;
    int lineStyle;
    SbBool lineWidthValid;
    int lineWidth;
    SbBool transparencyValid;
    float transparency;
    SbBool colorOverrideValid;
    SbBool colorOverride;
    SbBool colorValid;
    SbColor color;
    SbBool selectedColorValid;
    SbColor selectedColor;
    SbBool highlightedColorValid;
    SbColor highlightedColor;
    SbBool ghostedColorValid;
    SbColor ghostedColor;
};

struct BOBOL_EXPORT BObolRealizedShapeSummary {
    enum ShapeKind {
	SHAPE_UNKNOWN = 0,
	SHAPE_VLIST = 1,
	SHAPE_MESH = 2
    };

    BObolRealizedShapeSummary(void);

    SbBool valid;
    int shapeKind;
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    int ownerDrawMode;
    uint32_t ownerSourceRevision;
    uint32_t ownerInputsRevision;
    uint32_t ownerViewRevision;
    uint32_t ownerRealizedRevision;
    uint32_t ownerRealizedSourceRevision;
    uint32_t ownerRealizedInputsRevision;
    uint32_t ownerRealizedViewRevision;
    int ownerRealizationStatus;
    SbString ownerRealizationDiagnostic;
    SbString ownerRealizationIdentity;
    SbBool ownerSourceStale;
    uint32_t ownerStaleReason;
    SbString displayName;
    SbString geometryName;
    SbString cacheIdentity;
    SbString sourceIdentity;
    SbBool databaseIntent;
    SbBool overlayIntent;
    SbBool hudIntent;
    SbBool localSource;
    SbBool sharedSource;
    SbBool nonDatabaseSource;
    int drawMode;
    SbString recordRole;
    SbString geometryKind;
    int regionId;
    int airCode;
    int materialId;
    int los;
    SbBool materialColorValid;
    SbColor materialColor;
    SbString materialShader;
    uint32_t materialRevision;
    SbBool visible;
    SbBool selectable;
    SbBool selected;
    SbBool highlighted;
    SbBool ghosted;
    SbBool hiddenLine;
    SbBool editEmphasis;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbString editIntentId;
    SbString editIntentRole;
    uint32_t lodPolicy;
    SbBool colorOverride;
    SbColor color;
    SbBool lodAvailable;
    int lodActiveLevel;
    uint64_t lodFaceCount;
    uint64_t lodPointCount;
    uint64_t lodOriginalPointCount;
    uint64_t lodNormalCount;
    SbBool lodHasSnappedPoints;
    SbBool lodHasNormals;
    SbVec3f lodBoundsMin;
    SbVec3f lodBoundsMax;
    int pointCount;
    int commandCount;
    int segmentCount;
    int pointPrimitiveCount;
    int triangleCount;
    int indexCount;
    SbBool boundsValid;
    SbBox3f bounds;
};

/** Immutable geometry and semantic state for one compact CAD occurrence. */
struct BOBOL_EXPORT BObolCompactOccurrence {
    BObolCompactOccurrence(void);

    std::shared_ptr<const Obol::PartGeometry> geometry;
    BObolRealizedShapeSummary summary;
    /* Maps geometry coordinates into the occurrence's object-local frame. */
    SbMatrix geometryTransform;
    SbMatrix localTransform;
    SbBool lodBacked;
    SbBool sourceMeshRequestValid;
    BObolSourceMeshRequest sourceMeshRequest;
    uint32_t occurrenceIndex;
    int booleanOperation;
};

struct BOBOL_EXPORT BObolRealizedMaterialSummary {
    BObolRealizedMaterialSummary(void);

    SbBool valid;
    SbString sourcePath;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    int ownerDrawMode;
    uint32_t ownerSourceRevision;
    uint32_t ownerInputsRevision;
    uint32_t ownerViewRevision;
    uint32_t ownerRealizedRevision;
    uint32_t ownerRealizedSourceRevision;
    uint32_t ownerRealizedInputsRevision;
    uint32_t ownerRealizedViewRevision;
    int ownerRealizationStatus;
    SbString ownerRealizationDiagnostic;
    SbString ownerRealizationIdentity;
    SbBool ownerSourceStale;
    uint32_t ownerStaleReason;
    SbString materialName;
    SbString parentName;
    SbString materialSource;
    int propertyCount;
};

struct BOBOL_EXPORT BObolSceneTreeSummary {
    enum NodeKind {
	NODE_UNKNOWN = 0,
	NODE_GROUP = 1,
	NODE_DATABASE_SOURCE = 2,
	NODE_VLIST_SHAPE = 3,
	NODE_MESH_SHAPE = 4,
	NODE_MATERIAL_OBJECT = 5,
	NODE_OTHER = 6
    };

    BObolSceneTreeSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool isGroup;
    SbBool isShape;
    SbBool isDatabaseSource;
    SbBool isMaterialObject;
    SbBool hasParent;
    int drawTreeDepth;
    int childCount;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    SbString displayName;
    SbString geometryName;
};

struct BOBOL_EXPORT BObolSceneDisplaySummary {
    BObolSceneDisplaySummary(void);

    SbBool valid;
    int nodeKind;
    SbBool isDatabaseSource;
    SbBool hasDrawIntent;
    SbString intentPath;
    int intentDrawMode;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    double transparency;
    int drawMode;
    SbBool materialValid;
    SbColor materialColor;
    uint32_t materialRevision;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

struct BOBOL_EXPORT BObolSceneMaterialSummary {
    BObolSceneMaterialSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool materialValid;
    uint32_t materialRevision;
    SbColor materialColor;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

struct BOBOL_EXPORT BObolSceneBoundsSummary {
    BObolSceneBoundsSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool boundsValid;
    SbBox3f bounds;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString ownerSourceInstanceKey;
    SbString path;
};

class BOBOL_EXPORT SoBRLDatabaseSource : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLDatabaseSource);

public:
    enum DrawMode {
	WIREFRAME = 0,
	SHADED = 1
    };

    enum RepresentationMode {
	REPRESENTATION_DEFAULT = -1,
	REPRESENTATION_WIRE = 0,
	REPRESENTATION_SHADED_BOTS = 1,
	REPRESENTATION_SHADED = 2,
	REPRESENTATION_EVAL_WIRE = 3,
	REPRESENTATION_HIDDEN_LINE = 4,
	REPRESENTATION_EVAL_POINTS = 5
    };

    enum RealizationStatus {
	UNREALIZED = 0,
	REALIZED = 1,
	FAILED = 2
    };

    enum StaleReason {
	STALE_NONE = 0,
	STALE_SOURCE = 1,
	STALE_INPUTS = 2,
	STALE_VIEW = 4,
	STALE_DATABASE = 8,
	STALE_DRAW = 16,
	STALE_TESSELLATION = 32
    };

    enum RealizationRoleFlag {
	REALIZATION_ROLE_NONE = 0,
	REALIZATION_ROLE_CSG = 1,
	REALIZATION_ROLE_MESH = 2,
	REALIZATION_ROLE_EXTERNAL = 4
    };

    enum MaterialPolicy {
	MATERIAL_INHERIT = 0,
	MATERIAL_DATABASE = 1
    };

    enum BooleanOperation {
	BOOLEAN_UNION = 0,
	BOOLEAN_SUBTRACT = 1,
	BOOLEAN_INTERSECT = 2
    };

    SoSFString instanceKey;
    SoSFString path;
    SoSFString parentInstanceKey;
    SoSFUInt32 occurrenceIndex;
    SoSFInt32 booleanOperation;
    SoSFString displayName;
    SoSFString representationKey;
    SoSFInt32 representationMode;
    SoSFBool auxiliarySource;
    SoSFEnum drawMode;
    SoSFBool visible;
    SoSFBool selected;
    SoSFBool highlighted;
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFUInt32 materialRevision;
    SoSFEnum materialPolicy;
    SoSFBool databaseMetadataValid;
    SoSFInt32 databaseRegionId;
    SoSFInt32 databaseAirCode;
    SoSFInt32 databaseMaterialId;
    SoSFInt32 databaseLos;
    SoSFBool databaseMaterialColorValid;
    SoSFColor databaseMaterialColor;
    SoSFString databaseMaterialShader;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFColor selectedColor;
    SoSFColor highlightedColor;
    SoSFColor ghostedColor;
    /* Local-to-scene placement metadata for this source instance.  Primary
     * geometry remains source-local; render traversal and summaries apply
     * this transform explicitly instead of baking it into point/index arrays.
     */
    SoSFBool drawMatrixValid;
    SoSFMatrix drawMatrix;
    SoSFBool drawCenterValid;
    SoSFVec3f drawCenter;
    SoSFBool drawSizeValid;
    SoSFFloat drawSize;
    SoSFBool sourceBoundsValid;
    SoSFVec3f sourceBoundsMin;
    SoSFVec3f sourceBoundsMax;
    SoSFFloat tessellationAbsTol;
    SoSFFloat tessellationRelTol;
    SoSFFloat tessellationNormTol;
    SoSFUInt32 lodBotThreshold;
    SoSFUInt32 sourceRevision;
    SoSFUInt32 inputsRevision;
    SoSFUInt32 viewRevision;
    SoSFUInt32 realizedRevision;
    SoSFUInt32 realizedSourceRevision;
    SoSFUInt32 realizedInputsRevision;
    SoSFUInt32 realizedViewRevision;
    SoSFEnum realizationStatus;
    SoSFString realizationDiagnostic;
    SoSFString realizationIdentity;
    SoSFInt32 realizationRoleFlags;
    SoSFBool realizationViewDependent;
    SoSFBool realizationCsgLodEnabled;
    SoSFBool realizationMeshLodEnabled;
    SoSFFloat realizationViewScale;
    SoSFFloat realizationLodScale;
    SoSFInt32 realizationViewWidth;
    SoSFInt32 realizationViewHeight;
    SoSFUInt32 realizationBotThreshold;
    SoSFFloat realizationCurveScale;
    SoSFFloat realizationPointScale;
    SoSFBool stale;
    SoSFUInt32 staleReason;

    SoBRLDatabaseSource(void);
    static void initClass(void);

    /**
     * Set the BRL-CAD database used for realization.
     *
     * The caller owns the db_i pointer.  SoBRLDatabaseSource stores a borrowed
     * pointer and does not close or otherwise manage the database lifetime.
     * The database must remain valid while this source can be realized.
     */
    void setDatabase(struct db_i *dbip);
    struct db_i *getDatabase(void) const;
    void setMeshLod(struct BObolMeshLod *lod);
    struct BObolMeshLod *getMeshLod(void) const;
    void clearMeshLod(void);
    int setMeshLodBounds(const SbVec3f &bmin, const SbVec3f &bmax);
    SbBool getMeshLodBounds(SbVec3f &bmin, SbVec3f &bmax) const;
    void clearMeshLodBounds(void);
    void configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);
    void configureDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);
    void configureDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourcePath,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode,
	struct db_i *database,
	int mode,
	uint32_t revision);
    /* Create an unattached, ref-counted source with copied field state and an
     * immutable database snapshot containing the draw root and its recursive
     * dependencies.  The caller must unref the returned source, db_close the
     * returned database, and delete snapshotPathOut when it is nonempty. */
    SoBRLDatabaseSource *createDetachedRealizationSource(
	struct db_i **databaseOut, SbString *snapshotPathOut = NULL) const;
    /* Transfer a completed detached compact registry into this live source.
     * Call only on the live source's owner thread. */
    int adoptDetachedCompactRealization(SoBRLDatabaseSource *detached);
    int retargetDatabaseSource(const char *sourcePath,
	uint32_t revision);
    int retargetDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	uint32_t revision);
    int setRepresentationState(const char *sourceRepresentationKey,
	int sourceRepresentationMode);
    int setHierarchyState(const char *sourceParentInstanceKey,
	uint32_t sourceOccurrenceIndex,
	int sourceBooleanOperation);

    void markStale(void);
    void markStale(uint32_t reason);
    int setDrawModeState(int drawMode);
    int setDisplayNameState(const char *name);
    int setMaterialPolicyState(int materialPolicy);
    int setRealizationState(int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setRealizationRoleFlags(int roleFlags);
    int setRealizationViewPolicy(SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int setDisplayState(SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool selected,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    int setDatabaseMetadataState(SbBool metadataValid,
	int regionId,
	int airCode,
	int materialId,
	int los,
	SbBool metadataMaterialColorValid,
	const SbColor &metadataMaterialColor,
	const SbString &metadataMaterialShader);
    int refreshMaterialColorFromDatabase(uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    int applyDisplayPatch(const BObolDatabaseSourceDisplayPatch &patch);
    /* Update local-to-scene placement metadata.  This does not mutate primary
     * geometry coordinates; database realization and external publication
     * must keep drawable arrays in source-local space.
     */
    int setPlacementState(SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setSourceBoundsState(SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax);
    void clearSourceBounds(void);
    SbBool getSourceBounds(SbBox3f &bounds) const;
    SbBool getEffectiveSourceBounds(SbBox3f &bounds) const;
    SbBool needsRealization(void) const;
    SbBool realizePrototypeWireframe(void);
    SbBool realizeDatabaseWireframe(void);
    SbBool realizeDatabaseMesh(void);
    int clearRealizedGeometry(SbBool preserveAuxiliary = TRUE);
    int clearExternalPrimaryGeometry(void);
    /* Publish source-local primary geometry for externally supplied source
     * realizations.  These APIs preserve source placement separately via
     * setPlacementState/drawMatrix.
     */
    int publishExternalLineSet(const BObolExternalLineSet &lineSet);
    int publishExternalPointSet(const BObolExternalPointSet &pointSet);
    int publishExternalTriangleMesh(
	const BObolExternalTriangleMesh &triangleMesh);
    int publishExternalAnnotation(
	const BObolExternalAnnotation &annotation);
    int publishPrimitiveWireframe(
	struct rt_db_internal *intern,
	const struct bg_tess_tol *ttol = NULL,
	const struct bn_tol *tol = NULL);
    SoBRLVListShape *getRealizedShape(void) const;
    SoBRLVListShape *getRealizedShape(int index) const;
    int getRealizedShapeCount(void) const;
    SbBool hasRealizedWireGeometry(void) const;
    SoBRLVListShape *findAuxiliaryVListShape(const char *name) const;
    SoBRLDatabaseSource *findAuxiliarySource(const char *sourcePath) const;
    int setAuxiliaryLineSet(const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int setAuxiliarySourceLineSet(const char *sourcePath,
	const char *auxDisplayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int clearAuxiliaryShapes(void);
    SoBRLMeshShape *getRealizedMesh(void) const;
    SoBRLMeshShape *getRealizedMesh(int index) const;
    int getRealizedMeshCount(void) const;
    SbBool hasRealizedMeshGeometry(void) const;
    SoBRLMaterialObject *getRealizedMaterialObject(void) const;
    SoBRLMaterialObject *getRealizedMaterialObject(int index) const;
    int getRealizedMaterialObjectCount(void) const;
    int setCompactOccurrence(const BObolCompactOccurrence &occurrence);
    int setCompactOccurrenceRegistry(
	const std::vector<BObolCompactOccurrence> &occurrences);
    SbBool hasCompactInstanceIndex(void) const;
    SbBool isCompactOccurrenceRegistry(void) const;
    int getCompactInstanceCount(void) const;
    /* Unique immutable payloads shared by compact occurrences. */
    int getCompactPartCount(void) const;
    SbBool getCompactInstanceHandle(int index,
	BObolCompactInstanceHandle &handle) const;
    SbBool getCompactOccurrence(int index,
	BObolCompactOccurrence &occurrence) const;
    /* Flatten visible compact wire instances for non-raster consumers such as
     * vector export without constructing per-occurrence Coin shape nodes. */
    SbBool copyCompactWireGeometry(std::vector<SbVec3f> &points,
	std::vector<int32_t> &commands) const;
    /* Copy one compact occurrence into source coordinates for a transient
     * selection/edit presentation.  The compact index remains authoritative
     * and retains no Coin shape for this operation. */
    SbBool copyCompactInstanceEditGeometry(
	const BObolCompactInstanceHandle &handle,
	std::vector<SbVec3f> &points,
	std::vector<int32_t> &commands,
	BObolCompactInstanceSummary &summary) const;
    SbBool getCompactInstanceSummary(
	const BObolCompactInstanceHandle &handle,
	BObolCompactInstanceSummary &summary) const;
    SbBool isCompactInstanceHandleValid(
	const BObolCompactInstanceHandle &handle) const;
    int getCompactInstanceCountForPath(const char *path,
	SbBool includeDescendants = TRUE) const;
    SbBool getCompactInstanceBoundsForPath(const char *path,
	SbBool includeDescendants, SbBox3f &bounds) const;
    int setCompactInstanceDisplayStateForPath(const char *path,
	SbBool includeDescendants,
	int visibleValid, SbBool visible,
	int selectedValid, SbBool selected,
	int highlightedValid, SbBool highlighted);
    int setCompactInstanceSelectableForPath(const char *path,
	SbBool includeDescendants, SbBool selectable);
    int setCompactInstanceRegionIdForPath(const char *path,
	SbBool includeDescendants, int regionId);
    int setCompactInstanceRegionMetadataForPath(const char *path,
	SbBool includeDescendants, int regionId, int airCode,
	int materialId, int los);
    int setCompactInstanceMetadataForPath(const char *path,
	SbBool includeDescendants, int regionId, int airCode,
	int materialId, int los, SbBool materialColorValid,
	const SbColor &materialColor, const SbString &materialShader);
    int setCompactSubtractLineStyle(int lineStyle);
    int refreshCompactObjectGeometry(const char *objectPath,
	uint32_t sourceRevision = 0);
    int prepareCompiledAssembly(void);
    SoBRLCadAssembly *compactViewLodAssembly(
	const std::vector<const BObolViewLodState::CadPayload *> &payloads)
	const;
    SbBool hasCompiledAssembly(void) const;
    int getCompiledAssemblyPartCount(void) const;
    int getCompiledAssemblyInstanceCount(void) const;
    SbBool getSummary(BObolDatabaseSourceSummary &summary) const;
    int getRealizedShapeSummaryCount(void) const;
    SbBool getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const;
    int getRealizedMaterialSummaryCount(void) const;
    SbBool getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const;
    SbBool getRealizedMaterialProperty(int materialIndex, int propertyIndex,
	SbString &groupOut, SbString &nameOut, SbString &valueOut) const;
    int getRealizedTreeSummaryCount(void) const;
    SbBool getRealizedTreeSummary(int index,
	BObolSceneTreeSummary &summary) const;
    int getRealizedDisplaySummaryCount(void) const;
    SbBool getRealizedDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const;
    int getRealizedSceneMaterialSummaryCount(void) const;
    SbBool getRealizedSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const;
    int getRealizedBoundsSummaryCount(void) const;
    SbBool getRealizedBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const;

protected:
    virtual ~SoBRLDatabaseSource(void);
    void GLRender(SoGLRenderAction *action) override;
    void GLRenderBelowPath(SoGLRenderAction *action) override;
    void callback(SoCallbackAction *action) override;
    void getBoundingBox(SoGetBoundingBoxAction *action) override;
    void rayPick(SoRayPickAction *action) override;

private:
    friend class SoBRLCadRenderBatch;
    friend class SoBRLExportAction;
    friend class SoBRLMeasureAction;
    friend class SoBRLSnapAction;
    friend SbBool bobol_database_source_realize_wireframe_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend SbBool bobol_database_source_realize_mesh_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend int bobol_database_source_realize_wireframe_compact_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend int bobol_database_source_realize_mesh_compact_with_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);
    friend void bobol_database_source_seed_realization_cache(
	    SoBRLDatabaseSource *source,
	    BObolDatabaseSourceRealizationCache *cache);

    static void fieldSensorCB(void *data, SoSensor *sensor);
    void attachFieldSensors(void);
    void detachFieldSensors(void);
    void syncRealizedShapeOwnerState(void);
    void syncCompactInstanceDisplayState(void);
    void rebuildCompactInstanceDisplayState(SbBool syncSourceState);
    void syncCompactInstancePlacementState(void);
    void seedCompactRealizationCache(
	BObolDatabaseSourceRealizationCache *cache) const;
    void clearCompiledAssembly(void);
    void ensureCompiledAssemblyChild(void);
    void markCompiledAssemblyDirty(void);
    void markCadBatchDirty(void);
    uint64_t cadBatchRevisionGet(void) const;
    void clearCompactInstanceIndex(void);
    void discardCompactInstanceHistory(void);
    void installCompactInstanceIndex(
	struct BObolCompactInstanceIndex *index,
	SbBool occurrenceRegistry);
    int syncCompiledAssembly(void);
    uint64_t cadBatchStructureSignature(void) const;
    uint64_t cadBatchStyleSignature(void) const;
    uint64_t cadBatchSemanticSignature(void) const;
    int appendCadRenderBatch(struct BObolCadBatchBuildState *state,
	SbBool includeGeometry, SbBool includeSemantics);
    const struct BObolCompactInstanceEntry *findCompactInstanceEntry(
	const BObolCompactInstanceHandle &handle) const;
    int exportCompactInstances(SoBRLExportAction *action,
	const SbMatrix &parentToWorld);
    int measureCompactInstances(SoBRLMeasureAction *action,
	const SbMatrix &parentToWorld);
    int snapCompactInstances(SoBRLSnapAction *action,
	const SbMatrix &parentToWorld);

    struct Impl;
    std::unique_ptr<Impl> d;
};

/**
 * Refresh effective database material colors for a retained source set.
 * Combination state and shared path prefixes are resolved once per sweep;
 * callers should prefer this over invoking the single-source database lookup
 * repeatedly.
 */
BOBOL_EXPORT int
bobol_database_sources_refresh_material_colors(
    SoBRLDatabaseSource *const *sources,
    size_t sourceCount,
    uint32_t materialRevision,
    struct db_i *dbip);

#endif /* BOBOL_BDATABASESOURCE_H */
