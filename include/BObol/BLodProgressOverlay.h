/*        B L O D P R O G R E S S O V E R L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodProgressOverlay.h */

#ifndef BOBOL_BLODPROGRESSOVERLAY_H
#define BOBOL_BLODPROGRESSOVERLAY_H

#include "BObol/BDefines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/nodes/SoSeparator.h>

class SoHUDKit;
struct BObolLodProgressPresentationStatus;

/** Retained, toolkit-neutral status card for progressive LoD work. */
class BOBOL_EXPORT SoBRLLodProgressOverlay : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLLodProgressOverlay);

public:
    SoSFString title;
    SoSFString detail;
    SoSFColor color;
    SoSFBool visible;
    SoSFBool terminal;
    SoSFBool terminalReady;

    SoBRLLodProgressOverlay(void);
    static void initClass(void);

    void setStatus(const BObolLodProgressPresentationStatus &status);
    SoHUDKit *rebuildGeometry(void);
    SoHUDKit *getHUDKit(void) const;

protected:
    virtual ~SoBRLLodProgressOverlay(void);
};

#endif /* BOBOL_BLODPROGRESSOVERLAY_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
