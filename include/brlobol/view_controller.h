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
#include "brlobol/scene_controller.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SbViewportRegion.h>
#include <atomic>
#include <stddef.h>
#include <stdint.h>

class BRLObolLodService;
class SoCamera;
class SoNode;
class SoRenderManager;
class SoViewport;
struct bg_line_layer_builder;
struct db_i;
struct rt_view_info;

/**
 * Narrow application-facing controller for an Obol-backed BRL-CAD view.
 *
 * The controller records view services that GUI and command layers need while
 * leaving scene hierarchy, camera fields, visibility, and renderable geometry
 * in Obol.  It is the migration target for code that currently carries a
 * legacy view/display-manager pair.
 */
class BRLOBOL_EXPORT BRLObolViewController {
public:
    BRLObolViewController(void);
    explicit BRLObolViewController(SoNode *root, SoCamera *camera = NULL);
    ~BRLObolViewController(void);

    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;
    SoNode *getRenderRoot(void) const;

    void setCamera(SoCamera *camera);
    SoCamera *getCamera(void) const;

    void setViewportRegion(const SbViewportRegion &region);
    const SbViewportRegion &getViewportRegion(void) const;
    void setViewportSize(unsigned int width, unsigned int height);
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

    SoBRLSceneController *getSceneController(void);
    const SoBRLSceneController *getSceneController(void) const;
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
    int removeEditPreview(const char *previewId);

    int replaceDatabaseSource(const char *sourcePath,
	struct db_i *dbip,
	int drawMode = SoBRLDatabaseSource::WIREFRAME,
	uint32_t sourceRevision = 0);
    int removeDatabaseSource(const char *sourcePath);
    int clearDatabaseSources(void);
    SoBRLDatabaseSource *getDatabaseSource(int index) const;
    int getDatabaseSourceCount(void) const;

private:
    void syncRenderManager(void);
    void advanceLodViewRevision(void);
    void advanceLodPolicyRevision(void);
    void syncLodViewSignature(SbBool advanceOnChange = TRUE);
    static void lodResultReadyCB(BRLObolLodService *service, void *userData);

    SoBRLSceneController sceneController;
    SoViewport *viewport;
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
    unsigned int lastLodVisitedMeshCount;
    unsigned int lastLodSubmittedTaskCount;
    unsigned int lastLodSkippedMeshCount;
    size_t lastLodResultCount;
    unsigned int lastLodMatchedResultCount;
    unsigned int lastLodAppliedResultCount;
    unsigned int lastLodRejectedResultCount;
    unsigned int lastLodUnmatchedResultCount;
    SbString lastLodDiagnostics;
};

#endif /* BRLOBOL_VIEW_CONTROLLER_H */
