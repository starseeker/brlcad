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

namespace {

class BRLOBOLNullContextManager : public SoDB::ContextManager {
public:
    virtual void *createOffscreenContext(unsigned int UNUSED(width), unsigned int UNUSED(height))
    {
	return NULL;
    }

    virtual SbBool makeContextCurrent(void *UNUSED(context))
    {
	return FALSE;
    }

    virtual void restorePreviousContext(void *UNUSED(context))
    {
    }

    virtual void destroyContext(void *UNUSED(context))
    {
    }
};

}

void
brlobol_init(SoDB::ContextManager *contextManager)
{
    static bool initialized = false;
    if (initialized) {
	if (contextManager)
	    SoDB::setContextManager(contextManager);
	return;
    }

    static BRLOBOLNullContextManager nullContextManager;
    SoDB::init(contextManager ? contextManager : &nullContextManager);

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
}
