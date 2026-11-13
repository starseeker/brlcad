/*      D A T A B A S E _ S O U R C E _ C A C H E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

/** @file database_source_cache.cpp
 *
 * Shared realization-cache ownership and mutation.  This unit deliberately
 * has no source-tree realization policy; it only manages reusable geometry.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BMeshShape.h"
#include "BObol/BVListShape.h"
#include "cad_publication_private.h"
#include "database_source_private.h"
#include "database_source_realization.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

template <typename ShapeT>
static void
unref_cached_geometry_map(BObolRealizationCacheMap<ShapeT *> &cached)
{
    for (typename BObolRealizationCacheMap<ShapeT *>::iterator it = cached.begin();
	 it != cached.end(); ++it) {
	if (it->second)
	    it->second->unref();
    }
    cached.clear();
}

template <typename ShapeT>
static void
store_cached_geometry_map(BObolRealizationCacheMap<ShapeT *> &cached,
			  const std::string &key,
			  ShapeT *shape)
{
    if (!shape)
	return;

    /*
     * Realization seeding rebuilds the per-action cache from already realized
     * shared geometry.  Keep the internal shared-geometry key aligned with the
     * exact cache key used for storage, including view-LoD policy suffixes.
     * Instance shapes keep their user-facing geometry names separately.
     */
    shape->geometryName = key.c_str();
    shape->cacheIdentity =
	record_identity_with_revision(key.c_str(), shape->sourceId.getValue());

    bobol_cache_geometry_reference(cached, key, shape);
}

BObolDatabaseSourceRealizationCache::BObolDatabaseSourceRealizationCache(void)
{
}

BObolDatabaseSourceRealizationCache::~BObolDatabaseSourceRealizationCache(void)
{
    unref_cached_geometry_map(this->sharedWireGeometry);
    unref_cached_geometry_map(this->sharedMeshVListGeometry);
    unref_cached_geometry_map(this->sharedMeshGeometry);
}

void
BObolDatabaseSourceRealizationCache::clear(void) noexcept
{
    BObolDatabaseSourceRealizationCache retired;
    this->swap(retired);
}

void
BObolDatabaseSourceRealizationCache::swap(
    BObolDatabaseSourceRealizationCache &other) noexcept
{
    this->sharedWireGeometry.swap(other.sharedWireGeometry);
    this->sharedWireBounds.swap(other.sharedWireBounds);
    this->sharedMeshVListGeometry.swap(other.sharedMeshVListGeometry);
    this->sharedMeshGeometry.swap(other.sharedMeshGeometry);
    this->sharedWireCadGeometry.swap(other.sharedWireCadGeometry);
    this->sharedMeshVListCadGeometry.swap(other.sharedMeshVListCadGeometry);
    this->sharedMeshCadGeometry.swap(other.sharedMeshCadGeometry);
}

std::string
bobol_realization_cache_object_name(const std::string &path)
{
    const size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

template <typename Map, typename Match>
static void
extract_realization_cache_entries(Map &values, Map &retired, Match match)
{
    for (auto it = values.begin(); it != values.end();) {
	auto current = it++;
	if (match(current->first))
	    retired.insert(values.extract(current));
    }
}

void
BObolDatabaseSourceRealizationCache::eraseObject(const std::string &path)
{
    const std::string name = bobol_realization_cache_object_name(path);
    if (name.empty())
	return;
    const std::string prefix = name + ":";
    auto matches = [&name, &prefix](const std::string &key) {
	return key == name || key.compare(0, prefix.size(), prefix) == 0;
    };
    BObolDatabaseSourceRealizationCache retired;
    extract_realization_cache_entries(this->sharedWireGeometry,
	retired.sharedWireGeometry, matches);
    extract_realization_cache_entries(this->sharedWireBounds,
	retired.sharedWireBounds, matches);
    extract_realization_cache_entries(this->sharedMeshVListGeometry,
	retired.sharedMeshVListGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshGeometry,
	retired.sharedMeshGeometry, matches);
    extract_realization_cache_entries(this->sharedWireCadGeometry,
	retired.sharedWireCadGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshVListCadGeometry,
	retired.sharedMeshVListCadGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshCadGeometry,
	retired.sharedMeshCadGeometry, matches);
}

template <typename Map>
class PreparedRealizationCacheMapRename {
public:
    PreparedRealizationCacheMapRename(Map &values,
	const std::string &oldName, const std::string &newName) : target(values)
    {
	try {
	    auto stage = [&](const auto &entry) {
		const std::string nextKey = newName +
		    entry.first.substr(oldName.size());
		this->removedKeys.insert(entry.first);
		this->removedKeys.insert(nextKey);
		if constexpr (std::is_pointer_v<typename Map::mapped_type>) {
		    if (entry.second)
			entry.second->ref();
		    try {
			this->next.emplace(nextKey, entry.second);
		    } catch (...) {
			if (entry.second)
			    entry.second->unref();
			throw;
		    }
		} else {
		    this->next.emplace(nextKey, entry.second);
		}
	    };
	    auto exact = values.find(oldName);
	    if (exact != values.end())
		stage(*exact);
	    const std::string prefix = oldName + ":";
	    for (auto it = values.lower_bound(prefix); it != values.end() &&
		it->first.compare(0, prefix.size(), prefix) == 0; ++it)
		stage(*it);
	} catch (...) {
	    this->release(this->next);
	    throw;
	}
    }
    ~PreparedRealizationCacheMapRename()
    {
	this->release(this->next);
	this->release(this->retired);
    }
    void commit() noexcept
    {
	for (const auto &key : this->removedKeys) {
	    auto found = this->target.find(key);
	    if (found != this->target.end())
		this->retired.insert(this->target.extract(found));
	}
	for (auto it = this->next.begin(); it != this->next.end();) {
	    auto current = it++;
	    this->target.insert(this->next.extract(current));
	}
    }
private:
    static void release(Map &values) noexcept
    {
	if constexpr (std::is_pointer_v<typename Map::mapped_type>)
	    for (const auto &entry : values)
		if (entry.second)
		    entry.second->unref();
	values.clear();
    }
    Map &target;
    Map next;
    Map retired;
    std::set<std::string> removedKeys;
};

struct BObolDatabaseSourceRealizationCache::ObjectRename::Impl {
    Impl(BObolDatabaseSourceRealizationCache &cache,
	const std::string &oldName, const std::string &newName) :
	wire(cache.sharedWireGeometry, oldName, newName),
	bounds(cache.sharedWireBounds, oldName, newName),
	meshWire(cache.sharedMeshVListGeometry, oldName, newName),
	mesh(cache.sharedMeshGeometry, oldName, newName),
	wireCad(cache.sharedWireCadGeometry, oldName, newName),
	meshWireCad(cache.sharedMeshVListCadGeometry, oldName, newName),
	meshCad(cache.sharedMeshCadGeometry, oldName, newName)
    {
    }
    void commit() noexcept
    {
	wire.commit(); bounds.commit(); meshWire.commit(); mesh.commit();
	wireCad.commit(); meshWireCad.commit(); meshCad.commit();
    }
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedWireGeometry)> wire;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedWireBounds)> bounds;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedMeshVListGeometry)> meshWire;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedMeshGeometry)> mesh;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedWireCadGeometry)> wireCad;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedMeshVListCadGeometry)> meshWireCad;
    PreparedRealizationCacheMapRename<decltype(BObolDatabaseSourceRealizationCache::sharedMeshCadGeometry)> meshCad;
};

BObolDatabaseSourceRealizationCache::ObjectRename::ObjectRename(
    std::unique_ptr<Impl> prepared) : d(std::move(prepared))
{
}

BObolDatabaseSourceRealizationCache::ObjectRename::~ObjectRename() = default;

void
BObolDatabaseSourceRealizationCache::ObjectRename::commit() noexcept
{
    this->d->commit();
}

std::unique_ptr<BObolDatabaseSourceRealizationCache::ObjectRename>
BObolDatabaseSourceRealizationCache::prepareObjectRename(
    const std::string &oldPath, const std::string &newPath)
{
    const std::string oldName = bobol_realization_cache_object_name(oldPath);
    const std::string newName = bobol_realization_cache_object_name(newPath);
    return std::unique_ptr<ObjectRename>(new ObjectRename(
	std::make_unique<ObjectRename::Impl>(*this, oldName, newName)));
}

void
BObolDatabaseSourceRealizationCache::renameObject(
    const std::string &oldPath, const std::string &newPath)
{
    const std::string oldName = bobol_realization_cache_object_name(oldPath);
    const std::string newName = bobol_realization_cache_object_name(newPath);
    if (oldName.empty() || newName.empty() || oldName == newName)
	return;
    auto prepared = this->prepareObjectRename(oldName, newName);
    prepared->commit();
}

void
BObolDatabaseSourceRealizationCache::eraseViewVariants(void) noexcept
{
    static const char marker[] = ":view-lod:";
    auto matches = [](const std::string &key) {
	return key.find(marker) != std::string::npos;
    };
    BObolDatabaseSourceRealizationCache retired;
    extract_realization_cache_entries(this->sharedWireGeometry,
	retired.sharedWireGeometry, matches);
    extract_realization_cache_entries(this->sharedWireBounds,
	retired.sharedWireBounds, matches);
    extract_realization_cache_entries(this->sharedMeshVListGeometry,
	retired.sharedMeshVListGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshGeometry,
	retired.sharedMeshGeometry, matches);
    extract_realization_cache_entries(this->sharedWireCadGeometry,
	retired.sharedWireCadGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshVListCadGeometry,
	retired.sharedMeshVListCadGeometry, matches);
    extract_realization_cache_entries(this->sharedMeshCadGeometry,
	retired.sharedMeshCadGeometry, matches);
}

void
BObolDatabaseSourceRealizationCache::storeWireGeometry(
    const std::string &key,
    SoBRLVListShape *shape)
{
    store_cached_geometry_map(this->sharedWireGeometry, key, shape);
}

void
BObolDatabaseSourceRealizationCache::storeWireBounds(
    const std::string &key,
    const SbBox3f &bounds)
{
    if (key.empty() || bounds.isEmpty())
	return;

    this->sharedWireBounds[key] = bounds;
}

void
BObolDatabaseSourceRealizationCache::storeMeshVListGeometry(
    const std::string &key,
    SoBRLVListShape *shape)
{
    store_cached_geometry_map(this->sharedMeshVListGeometry, key, shape);
}

void
BObolDatabaseSourceRealizationCache::storeMeshGeometry(
    const std::string &key,
    SoBRLMeshShape *shape)
{
    store_cached_geometry_map(this->sharedMeshGeometry, key, shape);
}

static std::shared_ptr<const Obol::PartGeometry>
store_cached_part_geometry(
    BObolRealizationCacheMap<BObolCachedPartGeometry> &cache,
    const std::string &key, Obol::PartGeometryBuilder &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    if (key.empty())
	return std::shared_ptr<const Obol::PartGeometry>();

    const std::shared_ptr<const Obol::PartGeometry> sharedGeometry =
	bobol_cad_build_geometry(
	    std::move(geometry), "realization-cache insertion");
    if (!sharedGeometry)
	return std::shared_ptr<const Obol::PartGeometry>();

    BObolCachedPartGeometry &stored = cache[key];
    stored.geometry = sharedGeometry;
    stored.sourceType = sourceType ? sourceType : "";
    stored.geometryKind = geometryKind ? geometryKind : "";
    if (bounds)
	stored.bounds = *bounds;
    else
	stored.bounds = compact_part_geometry_bounds(stored.geometry);
    stored.geometryTransform = SbMatrix::identity();
    stored.viewDependentCsgGeometry = viewDependentCsgGeometry;
    stored.lodBacked = lodBacked;
    stored.sourceMeshRequestValid = sourceMeshRequest != NULL;
    if (sourceMeshRequest)
	stored.sourceMeshRequest = *sourceMeshRequest;
    return stored.geometry;
}

static std::shared_ptr<const Obol::PartGeometry>
store_cached_part_geometry_reference(
    BObolRealizationCacheMap<BObolCachedPartGeometry> &cache,
    const std::string &key, const std::shared_ptr<const Obol::PartGeometry> &geometry,
    const SbMatrix &geometryTransform, const char *sourceType,
    const char *geometryKind, const SbBox3f *bounds, bool lodBacked,
    const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    if (key.empty() || !geometry)
	return std::shared_ptr<const Obol::PartGeometry>();

    BObolCachedPartGeometry &stored = cache[key];
    stored.geometry = geometry;
    stored.sourceType = sourceType ? sourceType : "";
    stored.geometryKind = geometryKind ? geometryKind : "";
    if (bounds)
	stored.bounds = *bounds;
    else
	stored.bounds = database_source_transform_bounds(
	    compact_part_geometry_bounds(geometry), geometryTransform);
    stored.geometryTransform = geometryTransform;
    stored.viewDependentCsgGeometry = viewDependentCsgGeometry;
    stored.lodBacked = lodBacked;
    stored.sourceMeshRequestValid = sourceMeshRequest != NULL;
    if (sourceMeshRequest)
	stored.sourceMeshRequest = *sourceMeshRequest;
    return stored.geometry;
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeWireCadGeometry(
    const std::string &key, Obol::PartGeometryBuilder &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    return store_cached_part_geometry(this->sharedWireCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest, viewDependentCsgGeometry);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshVListCadGeometry(
    const std::string &key, Obol::PartGeometryBuilder &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    return store_cached_part_geometry(this->sharedMeshVListCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest, viewDependentCsgGeometry);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshCadGeometry(
    const std::string &key, Obol::PartGeometryBuilder &&geometry,
    const char *sourceType, const char *geometryKind, const SbBox3f *bounds,
    bool lodBacked, const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    return store_cached_part_geometry(this->sharedMeshCadGeometry, key,
	std::move(geometry), sourceType, geometryKind, bounds, lodBacked,
	sourceMeshRequest, viewDependentCsgGeometry);
}

std::shared_ptr<const Obol::PartGeometry>
BObolDatabaseSourceRealizationCache::storeMeshCadGeometryReference(
    const std::string &key,
    const std::shared_ptr<const Obol::PartGeometry> &geometry,
    const SbMatrix &geometryTransform, const char *sourceType,
    const char *geometryKind, const SbBox3f *bounds, bool lodBacked,
    const BObolSourceMeshRequest *sourceMeshRequest,
    bool viewDependentCsgGeometry)
{
    return store_cached_part_geometry_reference(this->sharedMeshCadGeometry,
	key, geometry, geometryTransform, sourceType, geometryKind, bounds,
	lodBacked, sourceMeshRequest, viewDependentCsgGeometry);
}

static const BObolCachedPartGeometry *
find_cached_part_geometry(
    const BObolRealizationCacheMap<BObolCachedPartGeometry> &cache,
    const std::string &key)
{
    auto found = cache.find(key);
    return found == cache.end() || !found->second.geometry ? NULL :
	&found->second;
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findWireCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedWireCadGeometry, key);
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findMeshVListCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedMeshVListCadGeometry, key);
}

const BObolCachedPartGeometry *
BObolDatabaseSourceRealizationCache::findMeshCadGeometry(
    const std::string &key) const
{
    return find_cached_part_geometry(this->sharedMeshCadGeometry, key);
}
