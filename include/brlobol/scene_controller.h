/*                S C E N E _ C O N T R O L L E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/scene_controller.h */

#ifndef BRLOBOL_SCENE_CONTROLLER_H
#define BRLOBOL_SCENE_CONTROLLER_H

#include "brlobol/defines.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbString.h>

class SoNode;

class BRLOBOL_EXPORT SoBRLSceneController {
public:
    SoBRLSceneController(void);
    /**
     * Create a controller over an application-owned Obol scene root.
     *
     * The controller retains the root with normal Obol reference counting, but
     * it does not own an authoritative hierarchy separate from Obol.
     */
    explicit SoBRLSceneController(SoNode *root);
    ~SoBRLSceneController(void);

    /**
     * Replace the retained scene root used by subsequent realization passes.
     */
    void setSceneRoot(SoNode *root);
    SoNode *getSceneRoot(void) const;

    SbBool realizePending(void);

    unsigned int getLastVisitedSourceCount(void) const;
    unsigned int getLastRealizedSourceCount(void) const;
    unsigned int getLastFailedSourceCount(void) const;
    const SbString &getLastDiagnostics(void) const;

private:
    SoNode *root;
    unsigned int lastVisitedSourceCount;
    unsigned int lastRealizedSourceCount;
    unsigned int lastFailedSourceCount;
    SbString lastDiagnostics;
};

#endif /* BRLOBOL_SCENE_CONTROLLER_H */
