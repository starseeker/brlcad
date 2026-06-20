/*                  R E A L I Z E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/realize_action.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

SO_ACTION_SOURCE(SoBRLRealizeAction);

SoBRLRealizeAction::SoBRLRealizeAction(void) :
    visitedSourceCount(0),
    realizedSourceCount(0),
    failedSourceCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLRealizeAction);
}

SoBRLRealizeAction::~SoBRLRealizeAction(void)
{
}

void
SoBRLRealizeAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLRealizeAction, SoAction);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLRealizeAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoGroup, SoBRLRealizeAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLDatabaseSource, SoBRLRealizeAction::databaseSourceAction);
}

unsigned int
SoBRLRealizeAction::getVisitedSourceCount(void) const
{
    return this->visitedSourceCount;
}

unsigned int
SoBRLRealizeAction::getRealizedSourceCount(void) const
{
    return this->realizedSourceCount;
}

unsigned int
SoBRLRealizeAction::getFailedSourceCount(void) const
{
    return this->failedSourceCount;
}

const SbString &
SoBRLRealizeAction::getDiagnostics(void) const
{
    return this->diagnostics;
}

void
SoBRLRealizeAction::beginTraversal(SoNode *node)
{
    this->visitedSourceCount = 0;
    this->realizedSourceCount = 0;
    this->failedSourceCount = 0;
    this->diagnostics = "";
    this->traverse(node);
}

void
SoBRLRealizeAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

void
SoBRLRealizeAction::appendDiagnostic(const SoBRLDatabaseSource *source)
{
    if (this->diagnostics.getLength() > 0)
	this->diagnostics += "\n";

    this->diagnostics += source ? source->path.getValue() : SbString("<unknown>");
    this->diagnostics += ": ";
    if (source && source->realizationDiagnostic.getValue().getLength() > 0) {
	this->diagnostics += source->realizationDiagnostic.getValue();
    } else {
	this->diagnostics += "realization failed";
    }
}

void
SoBRLRealizeAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLRealizeAction *realizeAction = static_cast<SoBRLRealizeAction *>(action);
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);

    realizeAction->visitedSourceCount++;
    if (source->needsRealization()) {
	if (source->getDatabase()) {
	    SbBool realized = FALSE;
	    if (source->drawMode.getValue() == SoBRLDatabaseSource::SHADED)
		realized = source->realizeDatabaseMesh();
	    else
		realized = source->realizeDatabaseWireframe();

	    if (realized)
		realizeAction->realizedSourceCount++;
	    else {
		source->realizationStatus = SoBRLDatabaseSource::FAILED;
		realizeAction->failedSourceCount++;
		realizeAction->appendDiagnostic(source);
	    }
	} else if (source->realizePrototypeWireframe()) {
	    realizeAction->realizedSourceCount++;
	} else {
	    source->realizationStatus = SoBRLDatabaseSource::FAILED;
	    realizeAction->failedSourceCount++;
	    realizeAction->appendDiagnostic(source);
	}
    }

    source->doAction(action);
}
