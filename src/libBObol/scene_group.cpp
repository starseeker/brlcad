/*                    S C E N E _ G R O U P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLodRealization.h"
#include "BObol/BSceneGroup.h"

SO_NODE_SOURCE(SoBRLSceneGroup);

SoBRLSceneGroup::SoBRLSceneGroup(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLSceneGroup);

    SO_NODE_ADD_FIELD(groupPath, (""));
    SO_NODE_ADD_FIELD(drawIntentValid, (FALSE));
    SO_NODE_ADD_FIELD(drawIntentPath, (""));
    SO_NODE_ADD_FIELD(drawMode, (BOBOL_LOD_DRAW_UNKNOWN));
    SO_NODE_ADD_FIELD(fallbackDrawMode, (BOBOL_LOD_DRAW_WIRE));
    SO_NODE_ADD_FIELD(overlayIntent, (FALSE));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selected, (FALSE));
    SO_NODE_ADD_FIELD(highlighted, (FALSE));
    SO_NODE_ADD_FIELD(lineStyle, (0));
    SO_NODE_ADD_FIELD(lineWidth, (0));
    SO_NODE_ADD_FIELD(transparency, (0.0f));
    SO_NODE_ADD_FIELD(colorOverride, (FALSE));
    SO_NODE_ADD_FIELD(color, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(materialColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialRevision, (0));
    SO_NODE_ADD_FIELD(revalidationRevision, (0));
}

SoBRLSceneGroup::~SoBRLSceneGroup(void)
{
}

void
SoBRLSceneGroup::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLSceneGroup, SoSeparator, "Separator");
}
