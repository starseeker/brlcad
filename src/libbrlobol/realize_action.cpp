/*                  R E A L I Z E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/realize_action.h"
#include "database_source_realization.h"
#include "performance_private.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

SO_ACTION_SOURCE(SoBRLRealizeAction);

SoBRLRealizeAction::SoBRLRealizeAction(void) :
    visitedSourceCount(0),
    realizedSourceCount(0),
    failedSourceCount(0),
    diagnostics(""),
    realizationCache(new BRLObolDatabaseSourceRealizationCache),
    seedingCache(FALSE),
    compactCadRealizationEnabled(FALSE)
{
    SO_ACTION_CONSTRUCTOR(SoBRLRealizeAction);
}

SoBRLRealizeAction::~SoBRLRealizeAction(void)
{
    delete this->realizationCache;
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
SoBRLRealizeAction::setCompactCadRealizationEnabled(SbBool enabled)
{
    this->compactCadRealizationEnabled = enabled ? TRUE : FALSE;
}

SbBool
SoBRLRealizeAction::getCompactCadRealizationEnabled(void) const
{
    return this->compactCadRealizationEnabled;
}

void
SoBRLRealizeAction::beginTraversal(SoNode *node)
{
    BRLObolPerformanceTimer totalTimer(BRLOBOL_PERF_REALIZE_TOTAL_US);
    if (totalTimer.active())
	brlobol_performance_counter_add(BRLOBOL_PERF_REALIZE_CALLS, 1);

    this->visitedSourceCount = 0;
    this->realizedSourceCount = 0;
    this->failedSourceCount = 0;
    this->diagnostics = "";
    if (this->realizationCache)
	this->realizationCache->clear();
    this->seedingCache = TRUE;
    int64_t phaseStart = brlobol_performance_time_now();
    this->traverse(node);
    if (phaseStart > 0) {
	const int64_t elapsed = brlobol_performance_time_now() - phaseStart;
	if (elapsed > 0)
	    brlobol_performance_counter_add(BRLOBOL_PERF_REALIZE_SEED_US,
		static_cast<uint64_t>(elapsed));
    }
    this->seedingCache = FALSE;
    phaseStart = brlobol_performance_time_now();
    this->traverse(node);
    if (phaseStart > 0) {
	const int64_t elapsed = brlobol_performance_time_now() - phaseStart;
	if (elapsed > 0)
	    brlobol_performance_counter_add(BRLOBOL_PERF_REALIZE_WALK_US,
		static_cast<uint64_t>(elapsed));
    }
    brlobol_performance_counter_add(BRLOBOL_PERF_SOURCES_VISITED,
	static_cast<uint64_t>(this->visitedSourceCount));
    brlobol_performance_counter_add(BRLOBOL_PERF_SOURCES_REALIZED,
	static_cast<uint64_t>(this->realizedSourceCount));
    brlobol_performance_counter_add(BRLOBOL_PERF_SOURCES_FAILED,
	static_cast<uint64_t>(this->failedSourceCount));
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

    if (realizeAction->seedingCache) {
	brlobol_database_source_seed_realization_cache(
	    source, realizeAction->realizationCache);
	source->doAction(action);
	return;
    }

    realizeAction->visitedSourceCount++;
    const int roleFlags = source->realizationRoleFlags.getValue();
    if (source->needsRealization() &&
	!(roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL)) {
	if (source->getDatabase()) {
	    SbBool realized = FALSE;
	    int compactRealized = 0;
	    if ((roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
		source->representationMode.getValue() ==
		SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ||
		source->representationMode.getValue() ==
		SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS ||
		source->drawMode.getValue() == SoBRLDatabaseSource::SHADED) {
		if (realizeAction->compactCadRealizationEnabled)
		    compactRealized =
			brlobol_database_source_realize_mesh_compact_with_cache(
			    source, realizeAction->realizationCache);
		if (compactRealized > 0)
		    realized = TRUE;
		else if (compactRealized == 0)
		    realized = brlobol_database_source_realize_mesh_with_cache(
				   source, realizeAction->realizationCache);
	    } else {
		if (realizeAction->compactCadRealizationEnabled)
		    compactRealized =
			brlobol_database_source_realize_wireframe_compact_with_cache(
			    source, realizeAction->realizationCache);
		if (compactRealized > 0)
		    realized = TRUE;
		else if (compactRealized == 0)
		    realized =
			brlobol_database_source_realize_wireframe_with_cache(
			    source, realizeAction->realizationCache);
	    }

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
