/*                    S C E N E _ G R O U P . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/scene_group.h */

#ifndef BRLOBOL_SCENE_GROUP_H
#define BRLOBOL_SCENE_GROUP_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

class BRLOBOL_EXPORT SoBRLSceneGroup : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLSceneGroup);

public:
    SoSFString groupPath;
    SoSFBool drawIntentValid;
    SoSFString drawIntentPath;
    SoSFInt32 drawMode;
    SoSFInt32 fallbackDrawMode;
    SoSFBool overlayIntent;
    SoSFBool visible;
    SoSFBool selected;
    SoSFBool highlighted;
    SoSFInt32 lineStyle;
    SoSFInt32 lineWidth;
    SoSFFloat transparency;
    SoSFBool colorOverride;
    SoSFColor color;
    SoSFBool materialColorValid;
    SoSFColor materialColor;
    SoSFUInt32 materialRevision;
    SoSFUInt32 revalidationRevision;

    SoBRLSceneGroup(void);
    static void initClass(void);

protected:
    virtual ~SoBRLSceneGroup(void);
};

#endif /* BRLOBOL_SCENE_GROUP_H */
