/*               H U D _ L A B E L _ O V E R L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/hud_label_overlay.h */

#ifndef BRLOBOL_HUD_LABEL_OVERLAY_H
#define BRLOBOL_HUD_LABEL_OVERLAY_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec2f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoHUDKit;
class SoHUDLabel;

class BRLOBOL_EXPORT SoBRLHUDLabelOverlay : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLHUDLabelOverlay);

public:
    SoSFString labelId;
    SoSFUInt32 sourceId;
    SoSFString text;
    SoSFVec2f position;
    SoSFColor color;
    SoSFFloat fontSize;
    SoSFBool visible;

    SoBRLHUDLabelOverlay(void);
    static void initClass(void);

    SoHUDLabel *rebuildGeometry(void);
    SoHUDKit *getHUDKit(void) const;
    SoHUDLabel *getHUDLabel(void) const;

protected:
    virtual ~SoBRLHUDLabelOverlay(void);

private:
    SoHUDLabel *hudLabel;
};

#endif /* BRLOBOL_HUD_LABEL_OVERLAY_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
