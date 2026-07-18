/*                  R E A L I Z E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BRealizeAction.h"
#include "database_source_realization.h"
#include "performance_private.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>

SO_ACTION_SOURCE(SoBRLRealizeAction);

struct BObolRealizationRepository::Residency {
    std::unordered_map<const SoBRLDatabaseSource *, std::set<std::string> >
	sourceObjects;
    std::unordered_map<std::string, size_t> objectReferences;
};

static std::set<std::string>
repository_source_objects(SoBRLDatabaseSource *source)
{
    std::set<std::string> objects;
    if (!source)
	return objects;

    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary summary;
	if (!source->getCompactInstanceHandle(i, handle) ||
	    !source->getCompactInstanceSummary(handle, summary))
	    continue;
	const char *name = summary.sourceName.getString();
	if (name && name[0])
	    objects.insert(name);
    }
    if (source->hasCompactInstanceIndex())
	return objects;
    for (int i = 0; i < source->getRealizedShapeSummaryCount(); i++) {
	BObolRealizedShapeSummary summary;
	if (!source->getRealizedShapeSummary(i, summary))
	    continue;
	const char *name = summary.sourceName.getString();
	if (name && name[0])
	    objects.insert(name);
    }
    return objects;
}

template <typename ResidencyType>
static void
repository_release_objects(ResidencyType *residency,
    BObolDatabaseSourceRealizationCache *cache,
    const std::set<std::string> &objects)
{
    if (!residency)
	return;
    for (const std::string &name : objects) {
	auto found = residency->objectReferences.find(name);
	if (found == residency->objectReferences.end())
	    continue;
	if (found->second > 1) {
	    found->second--;
	    continue;
	}
	residency->objectReferences.erase(found);
	if (cache)
	    cache->eraseObject(name);
    }
}

BObolRealizationRepository::BObolRealizationRepository(void) :
    cache(new BObolDatabaseSourceRealizationCache),
    residency(new Residency)
{
}

BObolRealizationRepository::~BObolRealizationRepository(void)
{
    delete this->cache;
    delete this->residency;
}

void
BObolRealizationRepository::clear(void)
{
    if (this->cache)
	this->cache->clear();
}

void
BObolRealizationRepository::invalidateObject(const char *name)
{
    if (this->cache && name && name[0])
	this->cache->eraseObject(name);
}

void
BObolRealizationRepository::renameObject(
    const char *oldName, const char *newName)
{
    if (!oldName || !oldName[0] || !newName || !newName[0])
	return;
    if (this->cache)
	this->cache->renameObject(oldName, newName);
    if (!this->residency || bu_strcmp(oldName, newName) == 0)
	return;

    auto oldCount = this->residency->objectReferences.find(oldName);
    if (oldCount != this->residency->objectReferences.end()) {
	this->residency->objectReferences[newName] += oldCount->second;
	this->residency->objectReferences.erase(oldCount);
    }
    for (auto &sourceEntry : this->residency->sourceObjects) {
	if (sourceEntry.second.erase(oldName))
	    sourceEntry.second.insert(newName);
    }
}

void
BObolRealizationRepository::invalidateViewVariants(void)
{
    if (this->cache)
	this->cache->eraseViewVariants();
}

void
BObolRealizationRepository::seedSource(SoBRLDatabaseSource *source)
{
    if (!source)
	return;
    const std::set<std::string> nextObjects = repository_source_objects(source);
    std::set<std::string> previousObjects;
    auto previous = this->residency->sourceObjects.find(source);
    if (previous != this->residency->sourceObjects.end())
	previousObjects = previous->second;

    std::set<std::string> removed;
    std::set_difference(previousObjects.begin(), previousObjects.end(),
	nextObjects.begin(), nextObjects.end(),
	std::inserter(removed, removed.end()));
    repository_release_objects(this->residency, this->cache, removed);
    for (const std::string &name : nextObjects) {
	if (!previousObjects.count(name))
	    this->residency->objectReferences[name]++;
    }
    this->residency->sourceObjects[source] = nextObjects;
    if (this->cache)
	bobol_database_source_seed_realization_cache(source, this->cache);
}

void
BObolRealizationRepository::releaseSource(SoBRLDatabaseSource *source)
{
    if (!this->residency || !source)
	return;
    auto found = this->residency->sourceObjects.find(source);
    if (found == this->residency->sourceObjects.end())
	return;
    repository_release_objects(this->residency, this->cache, found->second);
    this->residency->sourceObjects.erase(found);
}

SoBRLRealizeAction::SoBRLRealizeAction(void) :
    visitedSourceCount(0),
    realizedSourceCount(0),
    failedSourceCount(0),
    diagnostics(""),
    realizationCache(NULL),
    realizationRepository(new BObolRealizationRepository),
    ownsRealizationRepository(TRUE),
    seedingCache(FALSE),
    retainRealizationCache(FALSE)
{
    this->realizationCache = this->realizationRepository->cache;
    SO_ACTION_CONSTRUCTOR(SoBRLRealizeAction);
}

SoBRLRealizeAction::~SoBRLRealizeAction(void)
{
    if (this->ownsRealizationRepository)
	delete this->realizationRepository;
    this->realizationRepository = NULL;
    this->realizationCache = NULL;
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
SoBRLRealizeAction::setRetainRealizationCache(SbBool retain)
{
    this->retainRealizationCache = retain ? TRUE : FALSE;
    if (!this->retainRealizationCache)
	this->clearRealizationCache();
}

void
SoBRLRealizeAction::clearRealizationCache(void)
{
    if (this->realizationRepository)
	this->realizationRepository->clear();
}

void
SoBRLRealizeAction::invalidateRealizationObject(const char *name)
{
    if (this->realizationRepository)
	this->realizationRepository->invalidateObject(name);
}

void
SoBRLRealizeAction::setRealizationRepository(
    BObolRealizationRepository *repository)
{
    if (!repository || repository == this->realizationRepository)
	return;
    if (this->ownsRealizationRepository)
	delete this->realizationRepository;
    this->realizationRepository = repository;
    this->realizationCache = repository->cache;
    this->ownsRealizationRepository = FALSE;
}

void
SoBRLRealizeAction::beginTraversal(SoNode *node)
{
    BObolPerformanceTimer totalTimer(BOBOL_PERF_REALIZE_TOTAL_US);
    if (totalTimer.active())
	bobol_performance_counter_add(BOBOL_PERF_REALIZE_CALLS, 1);

    this->visitedSourceCount = 0;
    this->realizedSourceCount = 0;
    this->failedSourceCount = 0;
    this->diagnostics = "";
    if (this->realizationCache && !this->retainRealizationCache)
	this->realizationCache->clear();
    this->seedingCache = TRUE;
    int64_t phaseStart = bobol_performance_time_now();
    this->traverse(node);
    if (phaseStart > 0) {
	const int64_t elapsed = bobol_performance_time_now() - phaseStart;
	if (elapsed > 0)
	    bobol_performance_counter_add(BOBOL_PERF_REALIZE_SEED_US,
		static_cast<uint64_t>(elapsed));
    }
    this->seedingCache = FALSE;
    phaseStart = bobol_performance_time_now();
    this->traverse(node);
    if (phaseStart > 0) {
	const int64_t elapsed = bobol_performance_time_now() - phaseStart;
	if (elapsed > 0)
	    bobol_performance_counter_add(BOBOL_PERF_REALIZE_WALK_US,
		static_cast<uint64_t>(elapsed));
    }
    bobol_performance_counter_add(BOBOL_PERF_SOURCES_VISITED,
	static_cast<uint64_t>(this->visitedSourceCount));
    bobol_performance_counter_add(BOBOL_PERF_SOURCES_REALIZED,
	static_cast<uint64_t>(this->realizedSourceCount));
    bobol_performance_counter_add(BOBOL_PERF_SOURCES_FAILED,
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
	if (realizeAction->realizationRepository)
	    realizeAction->realizationRepository->seedSource(source);
	else
	    bobol_database_source_seed_realization_cache(
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
	    const int representation = source->representationMode.getValue();
	    const bool evaluated =
		representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE ||
		representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS;
	    if ((roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
		representation ==
		SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ||
		representation ==
		SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS ||
		source->drawMode.getValue() == SoBRLDatabaseSource::SHADED) {
		realized = evaluated ?
		    bobol_database_source_realize_mesh_with_cache(
			source, realizeAction->realizationCache) :
		    (bobol_database_source_realize_mesh_compact_with_cache(
			source, realizeAction->realizationCache) > 0 ? TRUE : FALSE);
	    } else {
		realized = evaluated ?
		    bobol_database_source_realize_wireframe_with_cache(
			source, realizeAction->realizationCache) :
		    (bobol_database_source_realize_wireframe_compact_with_cache(
			source, realizeAction->realizationCache) > 0 ? TRUE : FALSE);
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

    if (realizeAction->realizationRepository)
	realizeAction->realizationRepository->seedSource(source);

    source->doAction(action);
}
