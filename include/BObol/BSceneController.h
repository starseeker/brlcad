/*            B S C E N E C O N T R O L L E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BSceneController.h */

#ifndef BOBOL_BSCENECONTROLLER_H
#define BOBOL_BSCENECONTROLLER_H

#include "BObol/BDefines.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>

#include <stdint.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class SoNode;
class SoGroup;
class SoBRLDatabaseSource;
class BObolRealizationRepository;
struct BObolAuxiliaryLineSetDisplayState;
struct BObolDatabaseSourceDisplayPatch;
struct BObolExternalAnnotation;
struct BObolExternalLineSet;
struct BObolExternalPointSet;
struct BObolExternalTriangleMesh;
struct db_i;
struct BObolDatabaseSourceSummary;
struct BObolRealizedMaterialSummary;
struct BObolRealizedShapeSummary;
struct BObolSceneBoundsSummary;
struct BObolSceneDisplaySummary;
struct BObolSceneMaterialSummary;
struct BObolSceneTreeSummary;

struct BOBOL_EXPORT BObolDatabaseSourcePublishState {
    BObolDatabaseSourcePublishState(void);

    const char *sourceInstanceKey;
    const char *sourcePath;
    const char *sourceRepresentationKey;
    const char *targetGroupPath;
    struct db_i *database;
    int drawMode;
    int representationMode;
    SbBool sourceRevisionValid;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool colorOverride;
    SbColor color;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
    SbBool roleFlagsValid;
    int roleFlags;
    SbBool viewPolicyValid;
    SbBool viewDependent;
    SbBool csgLodEnabled;
    SbBool meshLodEnabled;
    float viewScale;
    float lodScale;
    int viewWidth;
    int viewHeight;
    uint32_t botThreshold;
    float curveScale;
    float pointScale;
    SbBool placementValid;
    SbBool drawMatrixValid;
    SbMatrix drawMatrix;
    SbBool drawCenterValid;
    SbVec3f drawCenter;
    SbBool drawSizeValid;
    float drawSize;
};

struct BOBOL_EXPORT BObolSceneSummary {
    BObolSceneSummary(void);

    SbBool valid;
    SbBool hasRoot;
    SbBool rootIsGroup;
    int rootChildCount;
    int databaseSourceCount;
    int nonDatabaseRootChildCount;
    uint64_t structuralRevision;
    uint64_t frameRevision;
    unsigned int lastVisitedSourceCount;
    unsigned int lastRealizedSourceCount;
    unsigned int lastFailedSourceCount;
    SbString lastDiagnostics;
};

class BOBOL_EXPORT SoBRLSceneController {
public:
    SoBRLSceneController(void);
    /**
     * Create a controller over an application-owned Obol scene root.
     *
     * The controller retains the root with normal Obol reference counting, but
     * it does not own an authoritative hierarchy separate from Obol.
     */
    explicit SoBRLSceneController(SoNode *root);
    ~SoBRLSceneController(void);

    /**
     * Replace the retained scene root used by subsequent realization passes.
     */
    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    uint64_t getStructuralRevision(void) const;
    uint64_t getFrameRevision(void) const;
    SbBool getSceneSummary(BObolSceneSummary &summary) const;

    SbBool realizePending(void);
    void shareRealizationRepository(SoBRLSceneController *source);
    void clearRealizationRepository(void);
    void invalidateRealizationViewVariants(void);
    void renameRealizationObject(const char *oldName, const char *newName);
    void beginSceneMutationBatch(size_t expectedDatabaseSources = 0,
	size_t expectedGroups = 0);
    void endSceneMutationBatch(void);

    SoGroup *findGroup(const char *groupPath) const;
    SoGroup *ensureGroup(const char *groupPath);
    int setGroupDrawIntent(const char *groupPath,
	const char *intentPath,
	int drawMode,
	int fallbackDrawMode,
	SbBool overlayIntent,
	uint32_t revalidationRevision);
    int setGroupDisplayState(const char *groupPath,
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
    int renameGroup(const char *groupPath, const char *newLeafName);
    int appendChildToGroup(const char *groupPath, SoNode *child);
    int removeChildFromGroup(const char *groupPath, SoNode *child);
    int eraseGroupSubpath(const char *parentGroupPath,
	const char *subpath);
    int removeGroup(const char *groupPath);
    int clearGroup(const char *groupPath);
    int getGroupChildCount(const char *groupPath) const;
    int getGroupDescendantGroupCount(const char *groupPath) const;
    int getGroupDatabaseSourceCount(const char *groupPath) const;
    SbBool getDatabaseSourceBounds(SbBox3f &bounds,
	SbBool padForAutoview) const;
    SbBool getSceneSubtreeBounds(const char *nodePath,
	SbBool includeOverlays,
	SbBox3f &bounds) const;

    SoNode *findShape(const char *shapePath) const;
    SoGroup *findShapeParent(const char *shapePath) const;
    int moveShapeToGroup(const char *shapePath, const char *groupPath);
    int removeShape(const char *shapePath);
    int setShapeDrawState(const char *shapePath,
	int drawMode,
	SbBool databaseIntent,
	SbBool overlayIntent,
	SbBool hudIntent);
    int setShapeDisplayState(const char *shapePath,
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
    int setShapeSourceState(const char *shapePath,
	const char *ownerSourcePath,
	uint32_t ownerSourceRevision,
	uint32_t ownerInputsRevision,
	uint32_t ownerViewRevision,
	uint32_t ownerRealizedRevision,
	uint32_t ownerRealizedSourceRevision,
	uint32_t ownerRealizedInputsRevision,
	uint32_t ownerRealizedViewRevision,
	int ownerRealizationStatus,
	const char *ownerRealizationDiagnostic,
	const char *ownerRealizationIdentity,
	SbBool ownerSourceStale,
	uint32_t ownerStaleReason);
    int setShapePlacementState(const char *shapePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int publishDatabaseSourceAuxiliaryLineSet(const char *sourcePath,
	const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceAuxiliarySourceLineSet(
	const char *sourcePath,
	const char *auxiliarySourcePath,
	const char *displayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceInstanceAuxiliaryLineSet(
	const char *sourceInstanceKey,
	const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    int publishDatabaseSourceInstanceAuxiliarySourceLineSet(
	const char *sourceInstanceKey,
	const char *auxiliarySourcePath,
	const char *displayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BObolAuxiliaryLineSetDisplayState *displayState = NULL);
    /* External primary geometry is source-local.  The database source's
     * placement state carries local-to-scene transforms separately.
     */
    int publishDatabaseSourceExternalLineSet(const char *sourcePath,
	const BObolExternalLineSet &lineSet);
    int publishDatabaseSourceInstanceExternalLineSet(
	const char *sourceInstanceKey,
	const BObolExternalLineSet &lineSet);
    int publishDatabaseSourceInstancePrimitiveWireframe(
	const char *sourceInstanceKey,
	struct rt_db_internal *intern,
	const struct bg_tess_tol *ttol = NULL,
	const struct bn_tol *tol = NULL);
    int publishDatabaseSourceExternalPointSet(const char *sourcePath,
	const BObolExternalPointSet &pointSet);
    int publishDatabaseSourceInstanceExternalPointSet(
	const char *sourceInstanceKey,
	const BObolExternalPointSet &pointSet);
    int publishDatabaseSourceExternalTriangleMesh(const char *sourcePath,
	const BObolExternalTriangleMesh &triangleMesh);
    int publishDatabaseSourceInstanceExternalTriangleMesh(
	const char *sourceInstanceKey,
	const BObolExternalTriangleMesh &triangleMesh);
    int publishDatabaseSourceExternalAnnotation(const char *sourcePath,
	const BObolExternalAnnotation &annotation);
    int publishDatabaseSourceInstanceExternalAnnotation(
	const char *sourceInstanceKey,
	const BObolExternalAnnotation &annotation);
    int clearDatabaseSourceExternalPrimaryGeometry(const char *sourcePath);
    int clearDatabaseSourceInstanceExternalPrimaryGeometry(
	const char *sourceInstanceKey);
    int clearDatabaseSourceAuxiliaryShapes(const char *sourcePath);
    int clearDatabaseSourceInstanceAuxiliaryShapes(
	const char *sourceInstanceKey);

    SoBRLDatabaseSource *getDatabaseSource(int index) const;
    int getDatabaseSourceCount(void) const;
    SoBRLDatabaseSource *findDatabaseSource(const char *sourcePath) const;
    SoBRLDatabaseSource *findDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    int replaceDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int replaceDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *sourcePath,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int replaceDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourcePath,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode,
	struct db_i *database,
	int drawMode,
	uint32_t sourceRevision);
    int publishDatabaseSourceInstance(
	const BObolDatabaseSourcePublishState &state);
    int renameDatabaseSource(const char *sourcePath,
	const char *newSourcePath,
	uint32_t sourceRevision);
    int renameDatabaseSourceInstance(const char *sourceInstanceKey,
	const char *newSourceInstanceKey,
	const char *newSourcePath,
	uint32_t sourceRevision);
    int setDatabaseSourceState(const char *sourcePath,
	SbBool sourceRevisionValid,
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
    int setDatabaseSourceInstanceState(const char *sourceInstanceKey,
	SbBool sourceRevisionValid,
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
    int setDatabaseSourceDisplayPatch(const char *sourcePath,
	const BObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceInstanceDisplayPatch(const char *sourceInstanceKey,
	const BObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceDisplayName(const char *sourcePath,
	const char *displayName);
    int setDatabaseSourceInstanceDisplayName(const char *sourceInstanceKey,
	const char *displayName);
    int setDatabaseSourceDrawMode(const char *sourcePath,
	int drawMode);
    int setDatabaseSourceInstanceDrawMode(const char *sourceInstanceKey,
	int drawMode);
    int setDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode);
    int setDatabaseSourceMaterialPolicy(const char *sourcePath,
	int materialPolicy);
    int setDatabaseSourceInstanceMaterialPolicy(const char *sourceInstanceKey,
	int materialPolicy);
    int refreshDatabaseSourceInstanceMaterialColorFromDatabase(
	const char *sourceInstanceKey,
	uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    int refreshDatabaseSourceMaterialColorsFromDatabase(
	uint32_t materialRevision,
	struct db_i *overrideDbip = NULL);
    int setDatabaseSourcePlacementState(const char *sourcePath,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setDatabaseSourceInstancePlacementState(const char *sourceInstanceKey,
	SbBool drawMatrixValid,
	const SbMatrix &drawMatrix,
	SbBool drawCenterValid,
	const SbVec3f &drawCenter,
	SbBool drawSizeValid,
	float drawSize);
    int setDatabaseSourceInstanceHierarchyState(const char *sourceInstanceKey,
	const char *parentInstanceKey,
	uint32_t occurrenceIndex,
	int booleanOperation);
    int setDatabaseSourceBoundsState(const char *sourcePath,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax);
    int setDatabaseSourceInstanceBoundsState(const char *sourceInstanceKey,
	SbBool boundsValid,
	const SbVec3f &boundsMin,
	const SbVec3f &boundsMax);
    int markDatabaseSourceStale(const char *sourcePath,
	uint32_t staleReason);
    int markDatabaseSourceInstanceStale(const char *sourceInstanceKey,
	uint32_t staleReason);
    int refreshDatabaseSourceInstanceObject(const char *sourceInstanceKey,
	const char *objectPath, uint32_t sourceRevision = 0);
    int setDatabaseSourceRealizationState(const char *sourcePath,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setDatabaseSourceInstanceRealizationState(
	const char *sourceInstanceKey,
	int realizationStatus,
	uint32_t realizedSourceRevision,
	uint32_t realizedInputsRevision,
	uint32_t staleReason,
	const char *diagnostic = NULL);
    int setDatabaseSourceRealizationRoleFlags(const char *sourcePath,
	int roleFlags);
    int setDatabaseSourceInstanceRealizationRoleFlags(
	const char *sourceInstanceKey,
	int roleFlags);
    int setDatabaseSourceRealizationViewPolicy(const char *sourcePath,
	SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int setDatabaseSourceInstanceRealizationViewPolicy(
	const char *sourceInstanceKey,
	SbBool viewDependent,
	SbBool csgLodEnabled,
	SbBool meshLodEnabled,
	float viewScale,
	float lodScale,
	int viewWidth,
	int viewHeight,
	uint32_t botThreshold,
	float curveScale,
	float pointScale);
    int moveDatabaseSourceToGroup(const char *sourcePath,
	const char *groupPath);
    int moveDatabaseSourceInstanceToGroup(const char *sourceInstanceKey,
	const char *groupPath);
    int removeDatabaseSource(const char *sourcePath);
    int removeDatabaseSourceInstance(const char *sourceInstanceKey);
    int clearDatabaseSources(void);
    SbBool getDatabaseSourceSummary(int index,
	BObolDatabaseSourceSummary &summary) const;
    SbBool getDatabaseSourceSummaryForPath(const char *sourcePath,
	BObolDatabaseSourceSummary &summary) const;
    int getDatabaseSourceInstanceCountForPath(const char *sourcePath) const;
    SbBool getDatabaseSourceInstanceSummaryForPath(const char *sourcePath,
	int instanceIndex, BObolDatabaseSourceSummary &summary) const;
    SbBool getDatabaseSourceSummaryForInstance(const char *sourceInstanceKey,
	BObolDatabaseSourceSummary &summary) const;
    int getRealizedShapeSummaryCount(void) const;
    SbBool getRealizedShapeSummary(int index,
	BObolRealizedShapeSummary &summary) const;
    int getRealizedMaterialSummaryCount(void) const;
    SbBool getRealizedMaterialSummary(int index,
	BObolRealizedMaterialSummary &summary) const;
    SbBool getRealizedMaterialProperty(int materialIndex, int propertyIndex,
	SbString &groupOut, SbString &nameOut, SbString &valueOut) const;
    int getSceneTreeSummaryCount(void) const;
    SbBool getSceneTreeSummary(int index,
	BObolSceneTreeSummary &summary) const;
    SbBool getSceneTreeSummaryForPath(const char *nodePath,
	BObolSceneTreeSummary &summary) const;
    SbBool getSceneChildTreeSummary(const char *nodePath,
	int childIndex,
	BObolSceneTreeSummary &summary) const;
    int getSceneDisplaySummaryCount(void) const;
    SbBool getSceneDisplaySummary(int index,
	BObolSceneDisplaySummary &summary) const;
    int getSceneMaterialSummaryCount(void) const;
    SbBool getSceneMaterialSummary(int index,
	BObolSceneMaterialSummary &summary) const;
    int getSceneBoundsSummaryCount(void) const;
    SbBool getSceneBoundsSummary(int index,
	BObolSceneBoundsSummary &summary) const;

    unsigned int getLastVisitedSourceCount(void) const;
    unsigned int getLastRealizedSourceCount(void) const;
    unsigned int getLastFailedSourceCount(void) const;
    const SbString &getLastDiagnostics(void) const;

private:
    void advanceFrameRevision(void);
    void advanceStructuralRevision(void);
    void clearDatabaseSourceIndex(void) const;
    void rebuildDatabaseSourceIndex(void) const;
    void indexSceneGroup(SoGroup *group) const;
    SoGroup *findIndexedGroup(const char *groupPath) const;
    void indexDatabaseSource(SoBRLDatabaseSource *source,
	SoGroup *parent = NULL) const;
    void unindexDatabaseSource(SoBRLDatabaseSource *source) const;
    SoBRLDatabaseSource *findIndexedDatabaseSource(
	const char *sourcePath) const;
    SbString databaseSourceInstanceKeyForPath(
	const char *sourcePath) const;
    SoBRLDatabaseSource *findIndexedDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    SoGroup *findIndexedDatabaseSourceInstanceParent(
	const char *sourceInstanceKey) const;
    SbBool databaseSourceSummaryForSource(SoBRLDatabaseSource *source,
	BObolDatabaseSourceSummary &summary) const;

    SoNode *root;
    uint64_t structuralRevision;
    uint64_t frameRevision;
    int mutationBatchDepth;
    SbBool mutationBatchStructuralRevisionPending;
    SbBool mutationBatchFrameRevisionPending;
    mutable SbBool databaseSourceIndexValid;
    mutable std::unordered_map<std::string, SoGroup *> groupPathIndex;
    mutable std::unordered_map<std::string, SoBRLDatabaseSource *> databaseSourcePathIndex;
    mutable std::unordered_map<std::string, std::vector<SoBRLDatabaseSource *> > databaseSourcePathInstancesIndex;
    mutable std::unordered_map<std::string, SoBRLDatabaseSource *> databaseSourceInstanceIndex;
    mutable std::unordered_map<std::string, SoGroup *> databaseSourceInstanceParentIndex;
    mutable std::vector<SoBRLDatabaseSource *> databaseSourceOrder;
    mutable std::unordered_map<SoBRLDatabaseSource *, size_t> databaseSourceOrderIndex;
    unsigned int lastVisitedSourceCount;
    unsigned int lastRealizedSourceCount;
    unsigned int lastFailedSourceCount;
    SbString lastDiagnostics;
    std::shared_ptr<BObolRealizationRepository> realizationRepository;
};

#endif /* BOBOL_BSCENECONTROLLER_H */
