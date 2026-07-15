/*                         I N I T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/init.h"

#include "brlobol/adc.h"
#include "brlobol/axes.h"
#include "brlobol/database_source.h"
#include "brlobol/edit_preview.h"
#include "brlobol/export_action.h"
#include "brlobol/grid.h"
#include "brlobol/headless_window_host.h"
#include "brlobol/hud_label_overlay.h"
#include "brlobol/image_plane.h"
#include "brlobol/image_source.h"
#include "brlobol/line_layer_overlay.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_update_action.h"
#include "brlobol/material_object.h"
#include "brlobol/measure_action.h"
#include "brlobol/mesh_lod_submit_action.h"
#include "brlobol/mesh_residency_action.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/pick_detail.h"
#include "brlobol/realize_action.h"
#include "brlobol/scene_group.h"
#include "brlobol/snap_action.h"
#include "brlobol/viewport_image.h"
#include "brlobol/view_lod.h"
#include "brlobol/vlist_shape.h"

#include <Inventor/SoDB.h>

void
brlobol_init(SoDB::ContextManager *contextManager)
{
    static bool initialized = false;
    if (initialized) {
	if (contextManager)
	    SoDB::setContextManager(contextManager);
	return;
    }
    if (!contextManager)
	contextManager = brlobol_headless_context_manager();
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
    (void)brlobol_headless_host_factory_register();
}

SoDB::ContextManager *
brlobol_headless_context_manager(void)
{
    static SoDB::ContextManager *manager = SoDB::createOSMesaContextManager();
    return manager;
}
