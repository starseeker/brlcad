/*                S C E N E _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/realize_action.h"
#include "brlobol/scene_controller.h"

#include <Inventor/nodes/SoNode.h>

SoBRLSceneController::SoBRLSceneController(void) :
    root(NULL),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
}

SoBRLSceneController::SoBRLSceneController(SoNode *sceneRoot) :
    root(NULL),
    lastVisitedSourceCount(0),
    lastRealizedSourceCount(0),
    lastFailedSourceCount(0),
    lastDiagnostics("")
{
    this->setSceneRoot(sceneRoot);
}

SoBRLSceneController::~SoBRLSceneController(void)
{
    this->setSceneRoot(NULL);
}

void
SoBRLSceneController::setSceneRoot(SoNode *sceneRoot)
{
    if (sceneRoot)
	sceneRoot->ref();
    if (this->root)
	this->root->unref();
    this->root = sceneRoot;
}

SoNode *
SoBRLSceneController::getSceneRoot(void) const
{
    return this->root;
}

SbBool
SoBRLSceneController::realizePending(void)
{
    this->lastVisitedSourceCount = 0;
    this->lastRealizedSourceCount = 0;
    this->lastFailedSourceCount = 0;
    this->lastDiagnostics = "";

    if (!this->root)
	return FALSE;

    SoBRLRealizeAction action;
    action.apply(this->root);
    this->lastVisitedSourceCount = action.getVisitedSourceCount();
    this->lastRealizedSourceCount = action.getRealizedSourceCount();
    this->lastFailedSourceCount = action.getFailedSourceCount();
    this->lastDiagnostics = action.getDiagnostics();
    return this->lastFailedSourceCount == 0;
}

unsigned int
SoBRLSceneController::getLastVisitedSourceCount(void) const
{
    return this->lastVisitedSourceCount;
}

unsigned int
SoBRLSceneController::getLastRealizedSourceCount(void) const
{
    return this->lastRealizedSourceCount;
}

unsigned int
SoBRLSceneController::getLastFailedSourceCount(void) const
{
    return this->lastFailedSourceCount;
}

const SbString &
SoBRLSceneController::getLastDiagnostics(void) const
{
    return this->lastDiagnostics;
}
