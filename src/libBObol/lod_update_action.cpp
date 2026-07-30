/*            L O D _ U P D A T E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BViewLod.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <stdlib.h>
#include <string.h>
#include <utility>

SO_ACTION_SOURCE(SoBRLLodUpdateAction);

SoBRLLodUpdateAction::SoBRLLodUpdateAction(void) :
    viewState(NULL),
    matchedResultCount(0),
    appliedResultCount(0),
    rejectedResultCount(0),
    unmatchedResultCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLLodUpdateAction);
}

SoBRLLodUpdateAction::~SoBRLLodUpdateAction(void)
{
}

void
SoBRLLodUpdateAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLLodUpdateAction, SoAction);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLLodUpdateAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoGroup, SoBRLLodUpdateAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLDatabaseSource,
			 SoBRLLodUpdateAction::databaseSourceAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape, SoBRLLodUpdateAction::meshShapeAction);
}

void
SoBRLLodUpdateAction::clearResults(void)
{
    this->results.clear();
    this->matched.clear();
}

void
SoBRLLodUpdateAction::addResult(const BObolLodResult &result)
{
    this->results.push_back(result);
}

void
SoBRLLodUpdateAction::addResult(BObolLodResult &&result)
{
    this->results.push_back(std::move(result));
}

void
SoBRLLodUpdateAction::setResults(
    const std::vector<BObolLodResult> &newResults)
{
    this->results = newResults;
    this->matched.clear();
}

size_t
SoBRLLodUpdateAction::drainService(BObolLodService &service,
				   size_t maxResults)
{
    std::vector<BObolLodResult> drained;
    size_t count = service.drainResults(drained, maxResults);

    this->results.insert(this->results.end(), drained.begin(), drained.end());
    this->matched.clear();
    return count;
}

void
SoBRLLodUpdateAction::setViewLodState(BObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

BObolViewLodState *
SoBRLLodUpdateAction::getViewLodState(void) const
{
    return this->viewState;
}

size_t
SoBRLLodUpdateAction::getResultCount(void) const
{
    return this->results.size();
}

unsigned int
SoBRLLodUpdateAction::getMatchedResultCount(void) const
{
    return this->matchedResultCount;
}

unsigned int
SoBRLLodUpdateAction::getAppliedResultCount(void) const
{
    return this->appliedResultCount;
}

unsigned int
SoBRLLodUpdateAction::getRejectedResultCount(void) const
{
    return this->rejectedResultCount;
}

unsigned int
SoBRLLodUpdateAction::getUnmatchedResultCount(void) const
{
    return this->unmatchedResultCount;
}

const SbString &
SoBRLLodUpdateAction::getDiagnostics(void) const
{
    return this->diagnostics;
}

static SbBool
lod_update_string_matches(const SbString &a, const SbString &b)
{
    return a.getLength() > 0 && b.getLength() > 0 &&
	   bu_strcmp(a.getString(), b.getString()) == 0 ? TRUE : FALSE;
}

static SbBool
lod_update_path_matches(const SbString &shapePath,
			const SbString &resultPath)
{
    if (lod_update_string_matches(shapePath, resultPath))
	return TRUE;

    if (shapePath.getLength() > 1 && resultPath.getLength() > 0 &&
	shapePath.getString()[0] == '/' &&
	bu_strcmp(shapePath.getString() + 1, resultPath.getString()) == 0)
	return TRUE;

    return FALSE;
}

static SbBool
lod_update_matches_shape(const SoBRLMeshShape *shape,
			 const BObolLodResult &result)
{
    if (!shape)
	return FALSE;

    if (lod_update_path_matches(shape->sourcePath.getValue(),
				result.request.objectPath))
	return TRUE;
    if (lod_update_string_matches(shape->sourceName.getValue(),
				  result.request.objectName))
	return TRUE;

    return FALSE;
}

static SbBool
lod_update_matches_source(const SoBRLDatabaseSource *source,
			  const BObolLodResult &result)
{
    if (!source)
	return FALSE;

    if (result.request.occurrenceKey.getLength() > 0) {
	if (!source->hasCompactInstanceIndex())
	    return FALSE;
	return source->hasCompactInstanceKey(
	    result.request.occurrenceKey.getString());
    }

    if (lod_update_path_matches(source->path.getValue(),
				result.request.objectPath))
	return TRUE;
    if (lod_update_string_matches(source->instanceKey.getValue(),
				  result.request.objectPath))
	return TRUE;
    if (source->hasCompactInstanceIndex() &&
	result.request.objectPath.getLength() > 0 &&
	source->getCompactInstanceCountForPath(
	    result.request.objectPath.getString(), FALSE) > 0)
	return TRUE;

    const char *path = source->path.getValue().getString();
    const char *leaf = path ? strrchr(path, '/') : NULL;
    leaf = (leaf && leaf[1]) ? leaf + 1 : path;
    if (leaf && result.request.objectName.getLength() > 0 &&
	bu_strcmp(leaf, result.request.objectName.getString()) == 0)
	return TRUE;

    return FALSE;
}

void
SoBRLLodUpdateAction::beginTraversal(SoNode *node)
{
    this->matched.assign(this->results.size(), FALSE);
    this->matchedResultCount = 0;
    this->appliedResultCount = 0;
    this->rejectedResultCount = 0;
    this->unmatchedResultCount = 0;
    this->diagnostics = "";

    this->traverse(node);
    this->finalizeUnmatchedDiagnostics();
}

void
SoBRLLodUpdateAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

void
SoBRLLodUpdateAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLLodUpdateAction *updateAction =
	static_cast<SoBRLLodUpdateAction *>(action);
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);

    if (!source->hasCompactInstanceIndex()) {
	source->doAction(action);
	return;
    }

    for (size_t i = 0; i < updateAction->results.size(); i++) {
	BObolLodResult &result = updateAction->results[i];
	if (updateAction->matched[i])
	    continue;
	if (!lod_update_matches_source(source, result))
	    continue;

	if (!updateAction->matched[i]) {
	    updateAction->matched[i] = TRUE;
	    updateAction->matchedResultCount++;
	}

	if (!updateAction->viewState) {
	    updateAction->rejectedResultCount++;
	    updateAction->appendDiagnostic(result,
					   "view-local compact CAD LoD update requires a view state");
	    continue;
	}

	if (updateAction->viewState->consumeSourceResult(source, result))
	    updateAction->appliedResultCount++;
	else {
	    updateAction->rejectedResultCount++;
	    if (getenv("BOBOL_LOD_TRACE_REJECTIONS")) {
		static unsigned int sourceRejectTraceCount = 0;
		if (sourceRejectTraceCount++ < 64)
		    bu_log("BObol LoD rejection reason=traversal-source "
			   "object=%s occurrence=%s level=%d requested=%d "
			   "status=%d kind=%d progressive=%d valid=%d "
			   "route=%llu source=%p diagnostic=%s\n",
			   result.request.objectName.getString(),
			   result.request.occurrenceKey.getString(),
			   result.geometry.activeLevel,
			   result.request.requestedLevel,
			   result.providerStatus, result.resultKind,
			   result.progressiveMesh ? 1 : 0,
			   result.progressiveMesh &&
				   result.progressiveMesh->isValid() ?
			       1 : 0,
			   static_cast<unsigned long long>(
			       result.request.sourceRoutingId),
			   static_cast<void *>(source),
			   result.diagnostic.getString());
	    }
	    updateAction->appendDiagnostic(result,
					   "view-local compact CAD LoD result rejected by source");
	}
    }

    source->doAction(action);
}

void
SoBRLLodUpdateAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLLodUpdateAction *updateAction =
	static_cast<SoBRLLodUpdateAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);

    for (size_t i = 0; i < updateAction->results.size(); i++) {
	BObolLodResult &result = updateAction->results[i];
	if (updateAction->matched[i])
	    continue;
	if (!lod_update_matches_shape(shape, result))
	    continue;

	if (!updateAction->matched[i]) {
	    updateAction->matched[i] = TRUE;
	    updateAction->matchedResultCount++;
	}

	if (!updateAction->viewState) {
	    updateAction->rejectedResultCount++;
	    updateAction->appendDiagnostic(result,
					   "view-local LoD update requires a view state");
	    continue;
	}

	if (updateAction->viewState->consumeDisplayResult(shape, result))
	    updateAction->appliedResultCount++;
	else {
	    updateAction->rejectedResultCount++;
	    updateAction->appendDiagnostic(result,
					   "view-local LoD result rejected by mesh");
	}
    }
}

void
SoBRLLodUpdateAction::appendDiagnostic(const BObolLodResult &result,
				       const char *message)
{
    if (this->diagnostics.getLength() > 0)
	this->diagnostics += "\n";

    SbString target = result.request.objectPath.getLength() > 0 ?
		      result.request.objectPath : result.request.objectName;
    this->diagnostics += target.getLength() > 0 ? target : SbString("<unknown>");
    this->diagnostics += ": ";
    this->diagnostics += message ? message : "LoD update diagnostic";
    if (result.diagnostic.getLength() > 0) {
	this->diagnostics += " (";
	this->diagnostics += result.diagnostic;
	this->diagnostics += ")";
    }
}

void
SoBRLLodUpdateAction::finalizeUnmatchedDiagnostics(void)
{
    for (size_t i = 0; i < this->matched.size(); i++) {
	if (this->matched[i])
	    continue;
	this->unmatchedResultCount++;
	this->appendDiagnostic(this->results[i],
			       "staged LoD result did not match any mesh");
    }
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
