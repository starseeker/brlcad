/*                         I N I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BInit.h"

#include "BObol/BADC.h"
#include "BObol/BAxes.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BEditPreview.h"
#include "BObol/BExportAction.h"
#include "BObol/BGrid.h"
#include "BObol/BHeadlessWindowHost.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BImagePlane.h"
#include "BObol/BImageSource.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodMeshShape.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMaterialObject.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshResidencyAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BPickDetail.h"
#include "BObol/BRealizeAction.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewportImage.h"
#include "BObol/BViewLod.h"
#include "BObol/BVListShape.h"

#include <Inventor/SoDB.h>

void
bobol_init(SoDB::ContextManager *contextManager)
{
    static bool initialized = false;
    if (initialized) {
	if (contextManager)
	    SoDB::setContextManager(contextManager);
	return;
    }
    if (!contextManager)
	contextManager = bobol_headless_context_manager();
    SoDB::init(contextManager);

    SoBRLPickDetail::initClass();
    SoBRLVListShape::initClass();
    SoBRLDatabaseSource::initClass();
    SoBRLEditPreview::initClass();
    SoBRLADC::initClass();
    SoBRLGrid::initClass();
    SoBRLAxes::initClass();
    SoBRLHUDLabelOverlay::initClass();
    SoBRLImageSource::initClass();
    SoBRLViewportImage::initClass();
    SoBRLImagePlane::initClass();
    SoBRLLineLayerOverlay::initClass();
    SoBRLMeshShape::initClass();
    SoBRLLodMeshShape::initClass();
    SoBRLMaterialObject::initClass();
    SoBRLSourceMeshPickAction::initClass();
    SoBRLRealizeAction::initClass();
    SoBRLViewLodElement::initClass();
    SoBRLViewLodGroup::initClass();
    SoBRLLodUpdateAction::initClass();
    SoBRLMeshLodSubmitAction::initClass();
    SoBRLMeshResidencyAction::initClass();
    SoBRLSnapAction::initClass();
    SoBRLMeasureAction::initClass();
    SoBRLExportAction::initClass();
    SoBRLSceneGroup::initClass();

    initialized = true;
    (void)bobol_headless_host_factory_register();
}

SoDB::ContextManager *
bobol_headless_context_manager(void)
{
    static SoDB::ContextManager *manager = SoDB::createOSMesaContextManager();
    return manager;
}
