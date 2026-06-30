/*                D A T A B A S E _ S O U R C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/database_source.h */

#ifndef BRLOBOL_DATABASE_SOURCE_H
#define BRLOBOL_DATABASE_SOURCE_H

#include "brlobol/defines.h"

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

class SoBRLVListShape;
class SoBRLMeshShape;
class SoBRLMaterialObject;
class SoFieldSensor;
class SoSensor;
struct db_i;

struct BRLOBOL_EXPORT BRLObolDatabaseSourceSummary {
    BRLObolDatabaseSourceSummary(void);

    SbBool valid;
    SbString path;
    SbBool hasParent;
    int drawTreeDepth;
    SbString parentGroupPath;
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
    float realizationViewScale;
    uint32_t realizationBotThreshold;
    float realizationCurveScale;
    float realizationPointScale;
    SbBool visible;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
    int materialPolicy;
    SbBool colorOverride;
    SbColor color;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
    SbBool stale;
    uint32_t staleReason;
    int realizedShapeCount;
    int realizedMeshCount;
    int realizedMaterialObjectCount;
};

struct BRLOBOL_EXPORT BRLObolRealizedShapeSummary {
    enum ShapeKind {
	SHAPE_UNKNOWN = 0,
	SHAPE_VLIST = 1,
	SHAPE_MESH = 2
    };

    BRLObolRealizedShapeSummary(void);

    SbBool valid;
    int shapeKind;
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
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
    SbString editIntentId;
    SbString editIntentRole;
    uint32_t lodPolicy;
    int pointCount;
    int commandCount;
    int segmentCount;
    int pointPrimitiveCount;
    int triangleCount;
    int indexCount;
    SbBool boundsValid;
    SbBox3f bounds;
};

struct BRLOBOL_EXPORT BRLObolRealizedMaterialSummary {
    BRLObolRealizedMaterialSummary(void);

    SbBool valid;
    SbString sourcePath;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    int ownerSourceIndex;
    SbString ownerSourcePath;
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

struct BRLOBOL_EXPORT BRLObolSceneTreeSummary {
    enum NodeKind {
	NODE_UNKNOWN = 0,
	NODE_GROUP = 1,
	NODE_DATABASE_SOURCE = 2,
	NODE_VLIST_SHAPE = 3,
	NODE_MESH_SHAPE = 4,
	NODE_MATERIAL_OBJECT = 5,
	NODE_OTHER = 6
    };

    BRLObolSceneTreeSummary(void);

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
    SbString path;
    SbString sourceName;
    SbString sourceType;
    uint32_t sourceId;
    SbString displayName;
    SbString geometryName;
};

struct BRLOBOL_EXPORT BRLObolSceneDisplaySummary {
    BRLObolSceneDisplaySummary(void);

    SbBool valid;
    int nodeKind;
    SbBool isDatabaseSource;
    SbBool hasDrawIntent;
    SbString intentPath;
    int intentDrawMode;
    SbBool visible;
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
    SbString path;
};

struct BRLOBOL_EXPORT BRLObolSceneMaterialSummary {
    BRLObolSceneMaterialSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool materialValid;
    uint32_t materialRevision;
    SbColor materialColor;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString path;
};

struct BRLOBOL_EXPORT BRLObolSceneBoundsSummary {
    BRLObolSceneBoundsSummary(void);

    SbBool valid;
    int nodeKind;
    SbBool boundsValid;
    SbBox3f bounds;
    int ownerSourceIndex;
    SbString ownerSourcePath;
    SbString path;
};

class BRLOBOL_EXPORT SoBRLDatabaseSource : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLDatabaseSource);

public:
    enum DrawMode {
	WIREFRAME = 0,
	SHADED = 1
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
	REALIZATION_ROLE_MESH = 2
    };

    enum MaterialPolicy {
	MATERIAL_INHERIT = 0,
	MATERIAL_DATABASE = 1
    };

    SoSFString path;
    SoSFEnum drawMode;
    SoSFBool visible;
    SoSFBool highlighted;
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFUInt32 materialRevision;
    SoSFEnum materialPolicy;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFBool drawMatrixValid;
    SoSFMatrix drawMatrix;
    SoSFBool drawCenterValid;
    SoSFVec3f drawCenter;
    SoSFBool drawSizeValid;
    SoSFFloat drawSize;
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
    SoSFFloat realizationViewScale;
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
    void configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision);
    int retargetDatabaseSource(const char *sourcePath,
	uint32_t revision);

    void markStale(void);
    void markStale(uint32_t reason);
    int setDrawModeState(int drawMode);
    int setMaterialPolicyState(int materialPolicy);
    int setRealizationState(int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setRealizationRoleFlags(int roleFlags);
    int setRealizationViewPolicy(SbBool viewDependent,
	float viewScale,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int setDisplayState(SbBool sourceRevisionValid,
	uint32_t sourceRevision,
	uint32_t inputsRevision,
	SbBool visible,
	SbBool highlighted,
	int lineStyle,
	int lineWidth,
	float transparency,
	SbBool colorOverride,
	const SbColor &color,
	SbBool materialColorValid,
	const SbColor &materialColor,
	uint32_t materialRevision);
    int setPlacementState(SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    SbBool needsRealization(void) const;
    SbBool realizePrototypeWireframe(void);
    SbBool realizeDatabaseWireframe(void);
    SbBool realizeDatabaseMesh(void);
    SoBRLVListShape *getRealizedShape(void) const;
    SoBRLVListShape *getRealizedShape(int index) const;
    int getRealizedShapeCount(void) const;
    SoBRLVListShape *findAuxiliaryVListShape(const char *name) const;
    int setAuxiliaryLineSet(const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count);
    int clearAuxiliaryShapes(void);
    SoBRLMeshShape *getRealizedMesh(void) const;
    SoBRLMeshShape *getRealizedMesh(int index) const;
    int getRealizedMeshCount(void) const;
    SoBRLMaterialObject *getRealizedMaterialObject(void) const;
    SoBRLMaterialObject *getRealizedMaterialObject(int index) const;
    int getRealizedMaterialObjectCount(void) const;
    SbBool getSummary(BRLObolDatabaseSourceSummary &summary) const;
    int getRealizedShapeSummaryCount(void) const;
    SbBool getRealizedShapeSummary(int index,
	BRLObolRealizedShapeSummary &summary) const;
    int getRealizedMaterialSummaryCount(void) const;
    SbBool getRealizedMaterialSummary(int index,
	BRLObolRealizedMaterialSummary &summary) const;
    SbBool getRealizedMaterialProperty(int materialIndex, int propertyIndex,
	SbString &groupOut, SbString &nameOut, SbString &valueOut) const;
    int getRealizedTreeSummaryCount(void) const;
    SbBool getRealizedTreeSummary(int index,
	BRLObolSceneTreeSummary &summary) const;
    int getRealizedDisplaySummaryCount(void) const;
    SbBool getRealizedDisplaySummary(int index,
	BRLObolSceneDisplaySummary &summary) const;
    int getRealizedSceneMaterialSummaryCount(void) const;
    SbBool getRealizedSceneMaterialSummary(int index,
	BRLObolSceneMaterialSummary &summary) const;
    int getRealizedBoundsSummaryCount(void) const;
    SbBool getRealizedBoundsSummary(int index,
	BRLObolSceneBoundsSummary &summary) const;

protected:
    virtual ~SoBRLDatabaseSource(void);

private:
    static void fieldSensorCB(void *data, SoSensor *sensor);
    void attachFieldSensors(void);
    void detachFieldSensors(void);
    void syncRealizedShapeOwnerState(void);

    struct db_i *dbip;
    SoFieldSensor *pathSensor;
    SoFieldSensor *drawModeSensor;
    SoFieldSensor *tessellationAbsTolSensor;
    SoFieldSensor *tessellationRelTolSensor;
    SoFieldSensor *tessellationNormTolSensor;
    SoFieldSensor *lodBotThresholdSensor;
    SoFieldSensor *sourceRevisionSensor;
    SoFieldSensor *inputsRevisionSensor;
    SoFieldSensor *viewRevisionSensor;
};

#endif /* BRLOBOL_DATABASE_SOURCE_H */
