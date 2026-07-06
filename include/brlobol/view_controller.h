/*                 V I E W _ C O N T R O L L E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/view_controller.h */

#ifndef BRLOBOL_VIEW_CONTROLLER_H
#define BRLOBOL_VIEW_CONTROLLER_H

#include "brlobol/defines.h"
#include "brlobol/database_source.h"
#include "brlobol/pick_detail.h"
#include "brlobol/scene_controller.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbRotation.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SbViewportRegion.h>
#include <atomic>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "vmath.h"

class BRLObolLodService;
class BRLObolViewLodState;
class BRLObolFeatureStore;
class BRLObolPolygonStore;
class BRLObolSelectionStore;
class SoBRLExportAction;
class SoBRLMeasureAction;
class SoBRLSnapAction;
class SoBRLViewLodGroup;
class SoCamera;
class SoGroup;
class SoNode;
class SoRenderManager;
class SoViewport;
struct bg_line_layer_builder;
struct db_i;
struct rt_view_info;

BRLOBOL_EXPORT SbMatrix brlobol_sbmatrix_from_brl_mat(const mat_t mat);
BRLOBOL_EXPORT SbRotation brlobol_camera_orientation_from_brl_rotation(
    const mat_t rotation);

/**
 * Narrow application-facing controller for an Obol-backed BRL-CAD view.
 *
 * The controller records view services that GUI and command layers need while
 * leaving scene hierarchy, camera fields, visibility, and renderable geometry
 * in Obol.  It is the migration target for code that currently carries a
 * legacy view/display-manager pair.
 */
class BRLOBOL_EXPORT BRLObolViewController
{
public:
    BRLObolViewController(void);
    explicit BRLObolViewController(SoNode *root, SoCamera *camera = NULL);
    ~BRLObolViewController(void);

    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    void setRenderSceneRoot(SoNode *root);
    SoNode *getRenderSceneRoot(void) const;
    SoNode *getRenderRoot(void) const;
    BRLObolViewLodState *getViewLodState(void) const;
    void clearViewLodState(void);

    void setCamera(SoCamera *camera);
    SoCamera *getCamera(void) const;

    void setViewportRegion(const SbViewportRegion &region);
    const SbViewportRegion &getViewportRegion(void) const;
    void setViewportSize(unsigned int width, unsigned int height);
    SbBool syncCameraFromRtViewContext(const void *viewCtx,
				       SbBool createCamera = TRUE);
    SbBool getRtViewInfo(struct rt_view_info *info) const;

    SbBool realizePending(void);

    unsigned int getLastVisitedSourceCount(void) const;
    unsigned int getLastRealizedSourceCount(void) const;
    unsigned int getLastFailedSourceCount(void) const;
    const SbString &getLastDiagnostics(void) const;

    void requestRender(const char *reason = NULL);
    void clearRenderRequest(void);
    SbBool consumeRenderRequest(SbString *reason = NULL);
    SbBool renderPending(SbBool clearWindow = TRUE,
			 SbBool clearZBuffer = TRUE,
			 SbString *reason = NULL);
    SbBool isRenderRequested(void) const;
    const SbString &getRenderReason(void) const;

    void setLodService(BRLObolLodService *service);
    BRLObolLodService *getLodService(void) const;
    void setLodAutoSubmit(SbBool enabled);
    SbBool isLodAutoSubmitEnabled(void) const;
    void setLodForcedLevel(int level);
    void clearLodForcedLevel(void);
    SbBool hasLodForcedLevel(void) const;
    int getLodForcedLevel(void) const;
    void setExactFullDetailBudget(uint64_t maxFaceCount,
				  uint64_t maxPointCount);
    uint64_t getMaxExactFullDetailFaceCount(void) const;
    uint64_t getMaxExactFullDetailPointCount(void) const;
    int consumeExportSourceFullDetail(SoBRLExportAction &exportAction,
				      uint64_t generation = 0,
				      int *submittedRequestCount = NULL);
    int consumeMeasureSourceFullDetail(SoBRLMeasureAction &measureAction,
				       uint64_t generation = 0,
				       int *submittedRequestCount = NULL);
    int consumeSnapSourceFullDetail(SoBRLSnapAction &snapAction,
				    uint64_t generation = 0,
				    int *submittedRequestCount = NULL);
    int prepareRtPickCaches(void);
    int getRtPickCacheCount(void) const;
    BRLObolRtPickCache *getRtPickCache(int index) const;
    uint32_t getRtPickCacheSourceRevision(int index) const;
    int pickSourceMeshExactRay(BRLObolSourceMeshPickResult &pick,
			       const SbVec3f &rayOrigin,
			       const SbVec3f &rayDirection,
			       uint64_t generation = 0,
			       int *submittedRequestCount = NULL);
    int pickRtExactRay(std::vector<BRLObolRtPickResult> &results,
		       const SbVec3f &rayOrigin,
		       const SbVec3f &rayDirection,
		       SbBool pickAll = FALSE);
    void clearRtPickCaches(void);
    void setMeshResidencyBudget(size_t maxResidentMeshBytes,
				SbBool evictDisplayPayloads = TRUE);
    void clearMeshResidencyBudget(void);
    SbBool hasMeshResidencyBudget(void) const;
    size_t getMaxResidentMeshBytes(void) const;
    SbBool isMeshResidencyDisplayEvictionEnabled(void) const;
    size_t evictMeshPayloadsToBudget(size_t maxBytes,
				     SbBool evictDisplayPayloads = TRUE);
    size_t getLastMeshBudgetInitialResidentBytes(void) const;
    size_t getLastMeshBudgetFinalResidentBytes(void) const;
    size_t getLastMeshBudgetFreedResidentBytes(void) const;
    size_t getLastMeshBudgetFreedFullDetailBytes(void) const;
    size_t getLastMeshBudgetFreedDisplayBytes(void) const;
    unsigned int getLastMeshBudgetVisitedMeshCount(void) const;
    unsigned int getLastMeshBudgetEvictedFullDetailMeshCount(void) const;
    unsigned int getLastMeshBudgetEvictedDisplayMeshCount(void) const;
    SbBool hasPendingLodResults(void) const;
    size_t processPendingLodResults(size_t maxResults = 0);
    int submitLodRequestsIfNeeded(SbBool refreshMissing = TRUE,
				  int reset = 0);
    int submitLodRequests(BRLObolLodService *service = NULL,
			  uint64_t generation = 0,
			  SbBool refreshMissing = TRUE,
			  int reset = 0);
    int applyLodResults(BRLObolLodService *service = NULL,
			size_t maxResults = 0);
    uint64_t getLodViewRevision(void) const;
    void setLodPolicyRevision(uint64_t revision);
    uint64_t getLodPolicyRevision(void) const;
    unsigned int getLastLodVisitedMeshCount(void) const;
    unsigned int getLastLodSubmittedTaskCount(void) const;
    unsigned int getLastLodSkippedMeshCount(void) const;
    size_t getLastLodResultCount(void) const;
    unsigned int getLastLodMatchedResultCount(void) const;
    unsigned int getLastLodAppliedResultCount(void) const;
    unsigned int getLastLodRejectedResultCount(void) const;
    unsigned int getLastLodUnmatchedResultCount(void) const;
    const SbString &getLastLodDiagnostics(void) const;
    size_t getActiveLodMeshPayloadCount(void) const;
    size_t getActiveLodProxyPayloadCount(int proxyKind) const;
    size_t getActiveLodCadPayloadCount(void) const;

    SoBRLSceneController *getSceneController(void);
    const SoBRLSceneController *getSceneController(void) const;
    BRLObolFeatureStore &features(void);
    const BRLObolFeatureStore &features(void) const;
    BRLObolPolygonStore &polygons(void);
    const BRLObolPolygonStore &polygons(void) const;
    BRLObolSelectionStore &selection(void);
    const BRLObolSelectionStore &selection(void) const;
    SoViewport *getViewport(void);
    const SoViewport *getViewport(void) const;
    SoRenderManager *getRenderManager(void);
    const SoRenderManager *getRenderManager(void) const;

    int replaceLineLayerOverlay(const char *overlayId,
				const struct bg_line_layer_builder *builder,
				uint32_t sourceId = 0,
				SbBool selectable = TRUE);
    int replaceHUDLabelOverlay(const char *labelId,
			       const char *text,
			       const SbVec2f &position,
			       const SbColor &color,
			       float fontSize = 12.0f,
			       uint32_t sourceId = 0);
    int removeHUDLabelOverlay(const char *labelId);
    int replaceEditPreview(const char *previewId,
			   const char *identity,
			   const SbVec3f *points,
			   const int32_t *commands,
			   int count,
			   uint32_t sourceRevision = 0,
			   uint32_t inputsRevision = 0);
    int replaceEditPreviewWithIntent(const char *previewId,
				     const char *identity,
				     const char *editIntentId,
				     const char *editIntentRole,
				     const SbVec3f *points,
				     const int32_t *commands,
				     int count,
				     uint32_t sourceRevision = 0,
				     uint32_t inputsRevision = 0);
    int removeEditPreview(const char *previewId);

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

    int replaceDatabaseSource(const char *sourcePath,
			      struct db_i *dbip,
			      int drawMode = SoBRLDatabaseSource::WIREFRAME,
			      uint32_t sourceRevision = 0);
    int replaceDatabaseSourceInstance(const char *sourceInstanceKey,
				      const char *sourcePath,
				      struct db_i *dbip,
				      int drawMode = SoBRLDatabaseSource::WIREFRAME,
				      uint32_t sourceRevision = 0);
    int setDatabaseSourceState(const char *sourcePath,
			       SbBool sourceRevisionValid,
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
    int setDatabaseSourceInstanceState(const char *sourceInstanceKey,
				       SbBool sourceRevisionValid,
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
    int setDatabaseSourceDisplayPatch(const char *sourcePath,
				      const BRLObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceInstanceDisplayPatch(const char *sourceInstanceKey,
	    const BRLObolDatabaseSourceDisplayPatch &patch);
    int setDatabaseSourceDisplayName(const char *sourcePath,
				     const char *displayName);
    int setDatabaseSourceInstanceDisplayName(const char *sourceInstanceKey,
	    const char *displayName);
    int setDatabaseSourceBoundsState(const char *sourcePath,
				     SbBool boundsValid,
				     const SbVec3f &boundsMin,
				     const SbVec3f &boundsMax);
    int setDatabaseSourceInstanceBoundsState(const char *sourceInstanceKey,
	    SbBool boundsValid,
	    const SbVec3f &boundsMin,
	    const SbVec3f &boundsMax);
    int setDatabaseSourceMaterialPolicy(const char *sourcePath,
					int materialPolicy);
    int setDatabaseSourceInstanceMaterialPolicy(const char *sourceInstanceKey,
	    int materialPolicy);
    int markDatabaseSourceStale(const char *sourcePath,
				uint32_t staleReason);
    int markDatabaseSourceInstanceStale(const char *sourceInstanceKey,
					uint32_t staleReason);
    int moveDatabaseSourceToGroup(const char *sourcePath,
				  const char *groupPath);
    int moveDatabaseSourceInstanceToGroup(const char *sourceInstanceKey,
					  const char *groupPath);
    int removeDatabaseSource(const char *sourcePath);
    int removeDatabaseSourceInstance(const char *sourceInstanceKey);
    int clearDatabaseSources(void);
    SoBRLDatabaseSource *getDatabaseSource(int index) const;
    int getDatabaseSourceCount(void) const;
    SoBRLDatabaseSource *findDatabaseSourceInstance(
	const char *sourceInstanceKey) const;
    SbBool getDatabaseSourceSummary(int index,
				    BRLObolDatabaseSourceSummary &summary) const;

private:
    void setViewportSceneGraphWithLod(SoNode *root);
    void syncRenderManager(void);
    void advanceLodViewRevision(void);
    void advanceLodPolicyRevision(void);
    void syncLodViewSignature(SbBool advanceOnChange = TRUE);
    size_t enforceMeshResidencyBudget(void);
    static void lodResultReadyCB(BRLObolLodService *service, void *userData);

    SoBRLSceneController sceneController;
    SoViewport *viewport;
    SoBRLViewLodGroup *renderLodRoot;
    BRLObolViewLodState *viewLodState;
    SoRenderManager *renderManager;
    SoCamera *activeCamera;
    SbViewportRegion viewportRegion;
    SbBool renderRequested;
    SbString renderReason;
    BRLObolLodService *lodService;
    uint64_t lodResultSubscriberId;
    std::atomic<int> lodResultsPending;
    SbBool lodAutoSubmit;
    uint64_t lodActiveGeneration;
    uint64_t lodLastSubmittedViewRevision;
    uint64_t lodLastSubmittedPolicyRevision;
    SbString lodLastSubmittedSourceSignature;
    SbString lodViewSignature;
    uint64_t lodViewRevision;
    uint64_t lodPolicyRevision;
    SbBool lodUseForcedLevel;
    int lodForcedLevel;
    uint64_t maxExactFullDetailFaceCount;
    uint64_t maxExactFullDetailPointCount;
    std::vector<BRLObolRtPickCache *> rtPickCaches;
    std::vector<SbString> rtPickCachePaths;
    std::vector<struct db_i *> rtPickCacheDatabases;
    std::vector<uint32_t> rtPickCacheSourceRevisions;
    SbBool meshResidencyBudgetEnabled;
    size_t maxResidentMeshBytes;
    SbBool meshResidencyEvictDisplayPayloads;
    size_t lastMeshBudgetInitialResidentBytes;
    size_t lastMeshBudgetFinalResidentBytes;
    size_t lastMeshBudgetFreedFullDetailBytes;
    size_t lastMeshBudgetFreedDisplayBytes;
    unsigned int lastMeshBudgetVisitedMeshCount;
    unsigned int lastMeshBudgetEvictedFullDetailMeshCount;
    unsigned int lastMeshBudgetEvictedDisplayMeshCount;
    unsigned int lastLodVisitedMeshCount;
    unsigned int lastLodSubmittedTaskCount;
    unsigned int lastLodSkippedMeshCount;
    size_t lastLodResultCount;
    unsigned int lastLodMatchedResultCount;
    unsigned int lastLodAppliedResultCount;
    unsigned int lastLodRejectedResultCount;
    unsigned int lastLodUnmatchedResultCount;
    SbString lastLodDiagnostics;
    BRLObolFeatureStore *featureStore;
    BRLObolPolygonStore *polygonStore;
    BRLObolSelectionStore *selectionStore;
};

#endif /* BRLOBOL_VIEW_CONTROLLER_H */
