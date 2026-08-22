/*          V I E W _ C O N T R O L L E R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H
#define LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H

#include "BObol/BDatabaseSource.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "view_lod_coordinator_state_private.h"
#include "bv.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>

#include <atomic>
#include <mutex>
#include <vector>

class SoCamera;
class SoClipPlane;
class SoGroup;
class SoNode;
class SoOffscreenRenderer;
class SoBRLViewLodGroup;

std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);
bool controller_lod_trace_enabled(const char *name, uint64_t viewRevision);
BObolLodPresentationPolicy::Population controller_lod_presentation_population(
    const BObolViewLodState *state, uint64_t sceneDomainRevision);
SbVec3f bobol_headlight_default_offset(void);
double controller_aspect_from_region(const SbViewportRegion &region);
void controller_configure_render_environment(SoViewport *viewport);
SoClipPlane *controller_clip_plane(SoViewport *viewport, SbBool minimum);

/* Uninstalled state shared by the responsibility-specific controller units.
 * It retains direct storage and one allocation: extracting scene forwarding
 * must not add virtual dispatch or accessors to LoD hot paths. */
struct BObolViewController::Impl : BObolLodCoordinator {
    explicit Impl(BObolViewController *owner) :
	headlightOffsetEye(bobol_headlight_default_offset()),
	featureStore(new BObolFeatureStore(owner)),
	polygonStore(new BObolPolygonStore(owner)),
	selectionStore(new BObolSelectionStore)
    {
    }

    BObolSceneController sceneController;
    SoViewport *viewport = new SoViewport;
    SoBRLViewLodGroup *renderLodRoot = NULL;
    SoNode *renderBatchRoot = NULL;
    SoNode *renderPresentationRoot = NULL;
    SoGroup *framebufferUnderlayRoot = NULL;
    SoGroup *framebufferInterlayRoot = NULL;
    SoGroup *framebufferOverlayRoot = NULL;
    BObolViewAttachment *viewAttachment = new BObolViewAttachment;
    SoRenderManager *renderManager = new SoRenderManager;
    SoOffscreenRenderer *imageRenderer = NULL;
    SoDB::ContextManager *imageRendererManager = NULL;
    SoCamera *activeCamera = NULL;
    SbViewportRegion viewportRegion = SbViewportRegion(1, 1);
    SbColor backgroundBottom = SbColor(0.0f, 0.0f, 0.0f);
    SbColor backgroundTop = SbColor(0.0f, 0.0f, 0.0f);
    BObolViewController::SoftwareWireMode softwareWireMode =
	BObolViewController::SOFTWARE_WIRE_AUTO;
    std::atomic<int> endpointGraphicalRenderingEnabled {1};
    SbBool transparencyEnabled = TRUE;
    SbBool antialiasingEnabled = FALSE;
    BObolViewController::LightingProfile lightingProfile =
	BObolViewController::LIGHTING_STUDIO;
    SbBool headlightEnabled = TRUE;
    SbBool headlightCameraTracked = TRUE;
    SbBool sceneLightsEnabled = FALSE;
    SbVec3f headlightOffsetEye;
    SbRotation lastCameraOrientation = SbRotation::identity();
    std::vector<BObolSceneLightRealization> sceneLights;
    double clipMinimum = BV_VIEW_MIN;
    double clipMaximum = BV_VIEW_MAX;
    mutable std::mutex renderRequestMutex;
    SbBool progressiveWorkPending = FALSE;
    SbBool renderRequested = FALSE;
    SbBool renderLodCapacityRelevant = FALSE;
    SbString renderReason = SbString("");
    uint64_t hostWorkRevision = 0;
    uint64_t renderRequestSerial = 0;
    std::mutex frameRequestMutex;
    BObolFrameRequestCallback frameRequestCallback = NULL;
    void *frameRequestUserData = NULL;
    mutable std::mutex presentationSyncMutex;
    BObolPresentationSyncCallback presentationSyncCallback = NULL;
    void *presentationSyncUserData = NULL;
    BObolFeatureStore *featureStore;
    BObolPolygonStore *polygonStore;
    BObolSelectionStore *selectionStore;
};

#endif /* LIBBOBOL_VIEW_CONTROLLER_PRIVATE_H */
