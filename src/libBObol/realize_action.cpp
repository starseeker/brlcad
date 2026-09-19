/*                  R E A L I Z E _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeshShape.h"
#include "BObol/BRealizeAction.h"
#include "BObol/BVListShape.h"
#include "database_source_realization.h"
#include "performance_private.h"

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/tools/SbModernUtils.h>

#include <algorithm>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

SO_ACTION_SOURCE(SoBRLRealizeAction);

struct BObolRealizationRepository::Residency {
    struct Source {
	std::set<std::string> objects;
	std::set<const BObolSceneController *> controllers;
    };
    std::unordered_map<const SoBRLDatabaseSource *, Source> sourceObjects;
    std::unordered_map<std::string, size_t> objectReferences;
    size_t controllerCount = 0;
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

BObolRealizationRepository::BObolRealizationRepository(void) :
    cache(new BObolDatabaseSourceRealizationCache),
    residency(new Residency)
{
}

BObolRealizationRepository::~BObolRealizationRepository(void)
{
    this->cache.reset();
    this->residency.reset();
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

struct BObolRealizationRepository::ObjectRename::Impl {
    Impl(BObolRealizationRepository &repository, const std::string &oldName,
	const std::string &newName) :
	cacheUpdate(repository.cache->prepareObjectRename(oldName, newName)),
	residency(*repository.residency), oldObject(oldName), newObject(newName)
    {
	size_t changedSources = 0;
	for (const auto &entry : this->residency.sourceObjects)
	    if (entry.second.objects.count(oldName))
		++changedSources;
	this->sources.reserve(changedSources);
	for (auto &entry : this->residency.sourceObjects) {
	    const bool ownsOld = entry.second.objects.count(oldName) != 0;
	    const bool ownsNew = entry.second.objects.count(newName) != 0;
	    if (ownsOld) {
		Source change{&entry.second, {}};
		if (!ownsNew) {
		    std::set<std::string> prepared;
		    auto inserted = prepared.insert(newName);
		    change.newName = prepared.extract(inserted.first);
		}
		this->sources.push_back(std::move(change));
	    }
	    if (ownsOld || ownsNew)
		++this->newReferences;
	}
	if (this->newReferences &&
	    !this->residency.objectReferences.count(this->newObject)) {
	    std::unordered_map<std::string, size_t> prepared;
	    auto inserted = prepared.emplace(this->newObject,
		this->newReferences);
	    this->newReference = prepared.extract(inserted.first);
	    this->residency.objectReferences.reserve(
		this->residency.objectReferences.size() + 1);
	}
    }
    void commit() noexcept
    {
	this->cacheUpdate->commit();
	for (auto &change : this->sources) {
	    change.source->objects.erase(this->oldObject);
	    if (!change.newName.empty())
		change.source->objects.insert(std::move(change.newName));
	}
	this->residency.objectReferences.erase(this->oldObject);
	auto current = this->residency.objectReferences.find(this->newObject);
	if (!this->newReferences) {
	    if (current != this->residency.objectReferences.end())
		this->residency.objectReferences.erase(current);
	} else if (current != this->residency.objectReferences.end()) {
	    current->second = this->newReferences;
	} else {
	    this->residency.objectReferences.insert(std::move(this->newReference));
	}
    }
    struct Source {
	BObolRealizationRepository::Residency::Source *source;
	std::set<std::string>::node_type newName;
    };
    std::unique_ptr<BObolDatabaseSourceRealizationCache::ObjectRename> cacheUpdate;
    BObolRealizationRepository::Residency &residency;
    std::string oldObject;
    std::string newObject;
    std::vector<Source> sources;
    size_t newReferences = 0;
    std::unordered_map<std::string, size_t>::node_type newReference;
};

BObolRealizationRepository::ObjectRename::ObjectRename(
    std::unique_ptr<Impl> prepared) : d(std::move(prepared))
{
}

BObolRealizationRepository::ObjectRename::~ObjectRename() = default;

void
BObolRealizationRepository::ObjectRename::commit() noexcept
{
    this->d->commit();
}

void
BObolRealizationRepository::renameObject(
    const char *oldName, const char *newName)
{
    if (!oldName || !oldName[0] || !newName || !newName[0])
	return;
    const std::string oldObject = bobol_realization_cache_object_name(oldName);
    const std::string newObject = bobol_realization_cache_object_name(newName);
    if (!this->cache || !this->residency || oldObject.empty() ||
	newObject.empty() || oldObject == newObject)
	return;
    auto update = this->prepareObjectRename(oldObject.c_str(),
	newObject.c_str());
    if (update)
	update->commit();
}

std::unique_ptr<BObolRealizationRepository::ObjectRename>
BObolRealizationRepository::prepareObjectRename(
    const char *oldName, const char *newName)
{
    if (!oldName || !oldName[0] || !newName || !newName[0] ||
	!this->cache || !this->residency)
	return nullptr;
    const std::string oldObject = bobol_realization_cache_object_name(oldName);
    const std::string newObject = bobol_realization_cache_object_name(newName);
    if (oldObject.empty() || newObject.empty() || oldObject == newObject)
	return nullptr;
    return std::unique_ptr<ObjectRename>(new ObjectRename(
	std::make_unique<ObjectRename::Impl>(*this, oldObject, newObject)));
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
    if (!source) return;
    auto update = this->prepareSourceSeed({source});
    update->commit();
}

void
BObolRealizationRepository::releaseSource(SoBRLDatabaseSource *source)
{
    if (!source || !this->residency->sourceObjects.count(source)) return;
    auto update = this->prepareSourceRelease({source});
    update->commit();
}

/* Transparent lookup avoids allocating a basename or variant prefix during
 * last-owner retirement. Exact names and their colon-delimited variants share
 * ownership; other objects with the same text prefix remain independent. */
template <typename Map, typename Visit>
static void
visit_cached_object(Map &values, std::string_view path, Visit visit)
{
    const auto slash = path.find_last_of('/');
    const auto name = slash == std::string_view::npos ? path : path.substr(slash + 1);
    if (name.empty()) return;
    for (auto it = values.lower_bound(name); it != values.end() && it->first.compare(0, name.size(), name) == 0;) {
	auto current = it++;
	if (current->first.size() == name.size() || current->first[name.size()] == ':') visit(current);
    }
}

/* Retain retired cache values through the enclosing publication's observers. */
template <typename Map>
class PreparedCacheRemoval {
public:
    PreparedCacheRemoval(Map &target, const std::set<std::string> &names) : values(target)
    {
	for (const auto &name : names)
	    visit_cached_object(target, name, [this](auto entry) { this->entries.push_back(entry); });
	std::sort(this->entries.begin(), this->entries.end(), [](auto a, auto b) { return a->first < b->first; });
	this->entries.erase(std::unique(this->entries.begin(), this->entries.end()), this->entries.end());
    }
    ~PreparedCacheRemoval()
    {
	if constexpr (std::is_pointer<typename Map::mapped_type>::value)
	    for (const auto &entry : this->retired) if (entry.second) entry.second->unref();
    }
    void commit() noexcept
    {
	for (auto it : this->entries) this->retired.insert(this->values.extract(it));
    }
private:
    Map &values;
    Map retired;
    std::vector<typename Map::iterator> entries;
};

template <typename Map>
static void
publish_cached_entries(Map &target, Map &candidate) noexcept
{
    for (auto it = candidate.begin(); it != candidate.end();) {
	auto next = it++;
	auto previous = target.extract(next->first);
	target.insert(candidate.extract(next));
	if (!previous.empty()) candidate.insert(std::move(previous));
    }
}

struct BObolRealizationRepository::SourceUpdate::Impl {
    struct Objects {
	Objects(Residency &state, const std::vector<SoBRLDatabaseSource *> &seeded,
	    const std::vector<SoBRLDatabaseSource *> &released, const std::vector<SourceMembership> &memberships) : residency(state)
	{
	    for (auto *source : seeded) if (source) this->stage(source).objects = repository_source_objects(source);
	    for (const auto &membership : memberships) {
		if (!membership.source || !membership.controller) continue;
		auto &next = this->stage(membership.source);
		if (membership.present) next.controllers.insert(membership.controller);
		else next.controllers.erase(membership.controller);
		this->membershipSources.insert(membership.source);
	    }
	    for (auto *source : released) if (source) {
		auto &next = this->stage(source);
		if (next.controllers.empty()) this->retired.insert(source);
	    }
	    for (auto *source : this->membershipSources)
		if (this->sourceEntries.at(source).controllers.empty()) this->retired.insert(source);
	    for (auto &entry : this->sourceEntries) {
		if (this->retired.count(entry.first)) entry.second.objects.clear();
		auto previous = state.sourceObjects.find(entry.first);
		if (previous != state.sourceObjects.end())
		    for (const auto &name : previous->second.objects) {
			auto &count = this->reference(name);
			if (count) --count;
		    }
		for (const auto &name : entry.second.objects) ++this->reference(name);
	    }
	    size_t added = 0;
	    for (const auto &entry : this->references) {
		if (!entry.second) this->names.insert(entry.first);
		else if (!state.objectReferences.count(entry.first)) ++added;
	    }
	    if (added) state.objectReferences.reserve(state.objectReferences.size() + added);
	    added = 0;
	    for (const auto &entry : this->sourceEntries)
		if (!this->retired.count(entry.first) && !state.sourceObjects.count(entry.first)) ++added;
	    if (added) state.sourceObjects.reserve(state.sourceObjects.size() + added);
	}
	Residency::Source &stage(SoBRLDatabaseSource *source)
	{
	    auto next = this->sourceEntries.try_emplace(source);
	    if (next.second) {
		this->owners.emplace_back(source);
		auto previous = this->residency.sourceObjects.find(source);
		if (previous != this->residency.sourceObjects.end()) next.first->second = previous->second;
	    }
	    return next.first->second;
	}
	size_t &reference(const std::string &name)
	{
	    auto found = this->residency.objectReferences.find(name);
	    return this->references.try_emplace(name, found == this->residency.objectReferences.end() ? 0 : found->second).first->second;
	}
	void commit() noexcept
	{
	    for (auto it = this->references.begin(); it != this->references.end();) {
		auto next = it++;
		auto previous = this->residency.objectReferences.find(next->first);
		if (!next->second) {
		    if (previous != this->residency.objectReferences.end()) this->residency.objectReferences.erase(previous);
		} else if (previous != this->residency.objectReferences.end()) previous->second = next->second;
		else this->residency.objectReferences.insert(this->references.extract(next));
	    }
	    for (auto it = this->sourceEntries.begin(); it != this->sourceEntries.end();) {
		auto next = it++;
		auto previous = this->residency.sourceObjects.find(next->first);
		if (this->retired.count(next->first)) {
		    if (previous != this->residency.sourceObjects.end()) this->residency.sourceObjects.erase(previous);
		} else if (previous != this->residency.sourceObjects.end()) {
		    previous->second.objects.swap(next->second.objects);
		    previous->second.controllers.swap(next->second.controllers);
		} else this->residency.sourceObjects.insert(this->sourceEntries.extract(next));
	    }
	}
	Residency &residency;
	decltype(Residency::sourceObjects) sourceEntries;
	decltype(Residency::objectReferences) references;
	std::unordered_set<const SoBRLDatabaseSource *> membershipSources, retired;
	std::vector<SbModernUtils::SoNodeRef> owners;
	std::set<std::string> names;
    } objects;

    Impl(BObolRealizationRepository &repository, const std::vector<SoBRLDatabaseSource *> &seeded,
	const std::vector<SoBRLDatabaseSource *> &released, const std::vector<SourceMembership> &memberships) :
	objects(*repository.residency, seeded, released, memberships), cache(*repository.cache),
	wire(cache.sharedWireGeometry, objects.names),
	bounds(cache.sharedWireBounds, objects.names),
	meshWire(cache.sharedMeshVListGeometry, objects.names),
	mesh(cache.sharedMeshGeometry, objects.names),
	wireCad(cache.sharedWireCadGeometry, objects.names),
	meshWireCad(cache.sharedMeshVListCadGeometry, objects.names),
	meshCad(cache.sharedMeshCadGeometry, objects.names)
    {
	for (auto *source : seeded)
	    if (source && !this->objects.retired.count(source)) bobol_database_source_seed_realization_cache(source, &this->candidate);
    }
    void commit() noexcept
    {
	wire.commit(); bounds.commit(); meshWire.commit(); mesh.commit();
	wireCad.commit(); meshWireCad.commit(); meshCad.commit(); objects.commit();
	publish_cached_entries(cache.sharedWireGeometry, candidate.sharedWireGeometry);
	publish_cached_entries(cache.sharedWireBounds, candidate.sharedWireBounds);
	publish_cached_entries(cache.sharedMeshVListGeometry, candidate.sharedMeshVListGeometry);
	publish_cached_entries(cache.sharedMeshGeometry, candidate.sharedMeshGeometry);
	publish_cached_entries(cache.sharedWireCadGeometry, candidate.sharedWireCadGeometry);
	publish_cached_entries(cache.sharedMeshVListCadGeometry, candidate.sharedMeshVListCadGeometry);
	publish_cached_entries(cache.sharedMeshCadGeometry, candidate.sharedMeshCadGeometry);
    }
    BObolDatabaseSourceRealizationCache &cache;
    BObolDatabaseSourceRealizationCache candidate;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedWireGeometry)> wire;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedWireBounds)> bounds;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedMeshVListGeometry)> meshWire;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedMeshGeometry)> mesh;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedWireCadGeometry)> wireCad;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedMeshVListCadGeometry)> meshWireCad;
    PreparedCacheRemoval<decltype(BObolDatabaseSourceRealizationCache::sharedMeshCadGeometry)> meshCad;
};

BObolRealizationRepository::SourceUpdate::SourceUpdate(std::unique_ptr<Impl> prepared) : d(std::move(prepared)) {}
BObolRealizationRepository::SourceUpdate::~SourceUpdate() = default;
void BObolRealizationRepository::SourceUpdate::commit() noexcept { this->d->commit(); }

std::unique_ptr<BObolRealizationRepository::SourceUpdate>
BObolRealizationRepository::prepareSourceRelease(const std::vector<SoBRLDatabaseSource *> &sources)
{
    return std::unique_ptr<SourceUpdate>(new SourceUpdate(std::make_unique<SourceUpdate::Impl>(*this, std::vector<SoBRLDatabaseSource *>(), sources, std::vector<SourceMembership>())));
}

std::unique_ptr<BObolRealizationRepository::SourceUpdate>
BObolRealizationRepository::prepareSourceSeed(const std::vector<SoBRLDatabaseSource *> &sources)
{
    return std::unique_ptr<SourceUpdate>(new SourceUpdate(std::make_unique<SourceUpdate::Impl>(*this, sources, std::vector<SoBRLDatabaseSource *>(), std::vector<SourceMembership>())));
}

std::unique_ptr<BObolRealizationRepository::SourceUpdate>
BObolRealizationRepository::prepareSourceMembership(const std::vector<SourceMembership> &memberships,
    const std::vector<SoBRLDatabaseSource *> &seeded)
{
    return std::unique_ptr<SourceUpdate>(new SourceUpdate(std::make_unique<SourceUpdate::Impl>(
	*this, seeded, std::vector<SoBRLDatabaseSource *>(), memberships)));
}

void
BObolRealizationRepository::attachController() noexcept
{
    ++this->residency->controllerCount;
}

void
BObolRealizationRepository::detachController(const BObolSceneController *controller) noexcept
{
    if (!this->residency->controllerCount) return;
    --this->residency->controllerCount;
    std::set<std::string> names;
    for (auto it = this->residency->sourceObjects.begin(); it != this->residency->sourceObjects.end();) {
	auto current = it++;
	auto &source = current->second;
	if (!source.controllers.erase(controller) || !source.controllers.empty()) continue;
	for (auto object = source.objects.begin(); object != source.objects.end();) {
	    auto next = object++;
	    auto count = this->residency->objectReferences.find(*next);
	    if (count != this->residency->objectReferences.end() && count->second && !--count->second) {
		this->residency->objectReferences.erase(count);
		names.insert(source.objects.extract(next));
	    }
	}
	this->residency->sourceObjects.erase(current);
    }
    BObolDatabaseSourceRealizationCache retired;
    auto release = [&names](auto &live, auto &previous) {
	for (const auto &name : names)
	    visit_cached_object(live, name, [&](auto entry) { previous.insert(live.extract(entry)); });
    };
    release(this->cache->sharedWireGeometry, retired.sharedWireGeometry);
    release(this->cache->sharedWireBounds, retired.sharedWireBounds);
    release(this->cache->sharedMeshVListGeometry, retired.sharedMeshVListGeometry);
    release(this->cache->sharedMeshGeometry, retired.sharedMeshGeometry);
    release(this->cache->sharedWireCadGeometry, retired.sharedWireCadGeometry);
    release(this->cache->sharedMeshVListCadGeometry, retired.sharedMeshVListCadGeometry);
    release(this->cache->sharedMeshCadGeometry, retired.sharedMeshCadGeometry);
}

bool
BObolRealizationRepository::hasSourceOwner(const BObolSceneController *controller, const SoBRLDatabaseSource *source) const noexcept
{
    auto found = this->residency->sourceObjects.find(source);
    return found != this->residency->sourceObjects.end() && found->second.controllers.count(controller);
}

bool
BObolRealizationRepository::acceptsSource(const SoBRLDatabaseSource *source) const noexcept
{
    if (!this->residency->controllerCount) return true;
    auto found = this->residency->sourceObjects.find(source);
    return found != this->residency->sourceObjects.end() && !found->second.controllers.empty();
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
    this->realizationCache = this->realizationRepository->cache.get();
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
    this->realizationCache = repository->cache.get();
    this->ownsRealizationRepository = FALSE;
}

void
SoBRLRealizeAction::stopSceneTraversal() noexcept
{
    /* Native child traversal keeps positional cursors. A scene edit retires
     * every enclosing cursor, including calls suspended by a nested apply. */
    for (auto *action = this; action; action = action->enclosingSceneAction)
	action->setTerminated(TRUE);
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
    if (!this->hasTerminated()) this->traverse(node);
    if (phaseStart > 0) {
	const int64_t elapsed = bobol_performance_time_now() - phaseStart;
	if (elapsed > 0)
	    bobol_performance_counter_add(BOBOL_PERF_REALIZE_SEED_US,
		static_cast<uint64_t>(elapsed));
    }
    this->seedingCache = FALSE;
    phaseStart = bobol_performance_time_now();
    if (!this->hasTerminated()) this->traverse(node);
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
    if (!action->hasTerminated() && node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

struct SoBRLRealizeAction::SourcePublication : BObolSourceRealizationEffects {
    explicit SourcePublication(SoBRLRealizeAction &owner) : action(owner) {}

    void prepare(const SoBRLDatabaseSource &source, bool success,
	const SbString &diagnostic, const std::vector<SoNode *> &children,
	const std::vector<SoNode *> &changedNodes) override
    {
	this->realized = success;
	if (!success) {
	    this->diagnostics = this->action.diagnostics;
	    if (this->diagnostics.getLength()) this->diagnostics += "\n";
	    this->diagnostics += source.path.getValue();
	    this->diagnostics += ": ";
	    this->diagnostics += diagnostic.getLength() ? diagnostic.getString() : "realization failed";
	}
	if (this->action.publicationEffects)
	    this->action.publicationEffects->prepare(source, success, diagnostic, children, changedNodes);
    }

    void commit(bool changed) noexcept override
    {
	if (this->realized) {
	    ++this->action.realizedSourceCount;
	} else {
	    ++this->action.failedSourceCount;
	    this->action.diagnostics = std::move(this->diagnostics);
	}
	if (this->action.publicationEffects)
	    this->action.publicationEffects->commit(changed);
    }

    void notify() override
    {
	if (this->action.publicationEffects)
	    this->action.publicationEffects->notify();
    }

    SoBRLRealizeAction &action;
    bool realized = false;
    SbString diagnostics;
};

void
SoBRLRealizeAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLRealizeAction *realizeAction = static_cast<SoBRLRealizeAction *>(action);
    if (realizeAction->hasTerminated()) return;
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
    SbModernUtils::SoNodeRef sourceOwner(source);

    if (realizeAction->seedingCache) {
	if (realizeAction->realizationRepository) {
	    if (realizeAction->realizationRepository->acceptsSource(source))
		realizeAction->realizationRepository->seedSource(source);
	} else
	    bobol_database_source_seed_realization_cache(
		source, realizeAction->realizationCache);
	if (!realizeAction->hasTerminated()) source->doAction(action);
	return;
    }

    realizeAction->visitedSourceCount++;
    const int roleFlags = source->realizationRoleFlags.getValue();
    if (source->needsRealization() &&
	!(roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL)) {
	SourcePublication publication(*realizeAction);
	if (source->getDatabase()) {
	    const int representation = source->representationMode.getValue();
	    const bool evaluated =
		representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE ||
		representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS;
	    const bool mesh = (roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
		representation == SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ||
		representation == SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS ||
		source->drawMode.getValue() == SoBRLDatabaseSource::SHADED;
	    if (evaluated) {
		const auto realize = mesh ? bobol_database_source_realize_mesh_with_cache :
		    bobol_database_source_realize_wireframe_with_cache;
		(void)realize(source, realizeAction->realizationCache, &publication);
	    } else {
		const auto realize = mesh ? bobol_database_source_realize_mesh_compact_with_cache :
		    bobol_database_source_realize_wireframe_compact_with_cache;
		(void)realize(source, realizeAction->realizationCache, nullptr, &publication);
	    }
	} else {
	    (void)bobol_database_source_realize_prototype(source, &publication);
	}
    }

    if (realizeAction->hasTerminated()) return;
    if (realizeAction->realizationRepository && realizeAction->realizationRepository->acceptsSource(source))
	realizeAction->realizationRepository->seedSource(source);

    if (!realizeAction->hasTerminated()) source->doAction(action);
}
