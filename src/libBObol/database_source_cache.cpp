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
#include <string>
#include <utility>
#include <vector>

template <typename ShapeT>
static void
unref_cached_geometry_map(std::map<std::string, ShapeT *> &cached)
{
    for (typename std::map<std::string, ShapeT *>::iterator it = cached.begin();
	 it != cached.end(); ++it) {
	if (it->second)
	    it->second->unref();
    }
    cached.clear();
}

template <typename ShapeT>
static void
store_cached_geometry_map(std::map<std::string, ShapeT *> &cached,
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

    typename std::map<std::string, ShapeT *>::iterator found =
	cached.find(key);
    if (found != cached.end()) {
	if (found->second == shape)
	    return;
	if (found->second)
	    found->second->unref();
    }

    shape->ref();
    cached[key] = shape;
}

BObolDatabaseSourceRealizationCache::BObolDatabaseSourceRealizationCache(void) :
    preserveCompactSourceOnFailure(false)
{
}

BObolDatabaseSourceRealizationCache::~BObolDatabaseSourceRealizationCache(void)
{
    this->clear();
}

void
BObolDatabaseSourceRealizationCache::clear(void)
{
    unref_cached_geometry_map(this->sharedWireGeometry);
    unref_cached_geometry_map(this->sharedMeshVListGeometry);
    unref_cached_geometry_map(this->sharedMeshGeometry);
    this->sharedWireBounds.clear();
    this->sharedWireCadGeometry.clear();
    this->sharedMeshVListCadGeometry.clear();
    this->sharedMeshCadGeometry.clear();
}

template <typename ValueT>
static void
erase_realization_cache_values(std::map<std::string, ValueT> &values,
	const std::string &name)
{
    values.erase(name);
    const std::string prefix = name + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0;)
	it = values.erase(it);
}

template <typename ShapeT>
static void
erase_realization_cache_shapes(std::map<std::string, ShapeT *> &values,
	const std::string &name)
{

    auto exact = values.find(name);
    if (exact != values.end()) {
	if (exact->second)
	    exact->second->unref();
	values.erase(exact);
    }
    const std::string prefix = name + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0;) {
	if (it->second)
	    it->second->unref();
	it = values.erase(it);
    }
}

void
BObolDatabaseSourceRealizationCache::eraseObject(const std::string &path)
{
    const size_t slash = path.find_last_of('/');
    const std::string name = slash == std::string::npos ? path :
	path.substr(slash + 1);
    if (name.empty())
	return;
    erase_realization_cache_shapes(this->sharedWireGeometry, name);
    erase_realization_cache_values(this->sharedWireBounds, name);
    erase_realization_cache_shapes(this->sharedMeshVListGeometry, name);
    erase_realization_cache_shapes(this->sharedMeshGeometry, name);
    erase_realization_cache_values(this->sharedWireCadGeometry, name);
    erase_realization_cache_values(this->sharedMeshVListCadGeometry, name);
    erase_realization_cache_values(this->sharedMeshCadGeometry, name);
}

static std::string
realization_cache_object_name(const std::string &path)
{
    const size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

template <typename ValueT>
static void
rename_realization_cache_values(std::map<std::string, ValueT> &values,
	const std::string &oldName, const std::string &newName)
{
    std::vector<std::string> keys;
    if (values.find(oldName) != values.end())
	keys.push_back(oldName);
    const std::string prefix = oldName + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0; ++it)
	keys.push_back(it->first);
    for (const std::string &key : keys) {
	auto it = values.find(key);
	if (it == values.end())
	    continue;
	const std::string nextKey = newName + key.substr(oldName.size());
	ValueT value = std::move(it->second);
	values.erase(it);
	values.erase(nextKey);
	values.emplace(nextKey, std::move(value));
    }
}

template <typename ShapeT>
static void
rename_realization_cache_shapes(std::map<std::string, ShapeT *> &values,
	const std::string &oldName, const std::string &newName)
{
    std::vector<std::string> keys;
    if (values.find(oldName) != values.end())
	keys.push_back(oldName);
    const std::string prefix = oldName + ":";
    for (auto it = values.lower_bound(prefix); it != values.end() &&
	it->first.compare(0, prefix.size(), prefix) == 0; ++it)
	keys.push_back(it->first);
    for (const std::string &key : keys) {
	auto it = values.find(key);
	if (it == values.end())
	    continue;
	const std::string nextKey = newName + key.substr(oldName.size());
	ShapeT *shape = it->second;
	values.erase(it);
	auto destination = values.find(nextKey);
	if (destination != values.end()) {
	    if (destination->second)
		destination->second->unref();
	    values.erase(destination);
	}
	values.emplace(nextKey, shape);
    }
}

void
BObolDatabaseSourceRealizationCache::renameObject(
    const std::string &oldPath, const std::string &newPath)
{
    const std::string oldName = realization_cache_object_name(oldPath);
    const std::string newName = realization_cache_object_name(newPath);
    if (oldName.empty() || newName.empty() || oldName == newName)
	return;
    rename_realization_cache_shapes(this->sharedWireGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedWireBounds, oldName, newName);
    rename_realization_cache_shapes(this->sharedMeshVListGeometry, oldName,
	newName);
    rename_realization_cache_shapes(this->sharedMeshGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedWireCadGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedMeshVListCadGeometry, oldName,
	newName);
    rename_realization_cache_values(this->sharedMeshCadGeometry, oldName,
	newName);
}

template <typename ValueT>
static void
erase_realization_cache_variant_values(
    std::map<std::string, ValueT> &values, const char *marker)
{
    for (auto it = values.begin(); it != values.end();) {
	if (it->first.find(marker) != std::string::npos)
	    it = values.erase(it);
	else
	    ++it;
    }
}

template <typename ShapeT>
static void
erase_realization_cache_variant_shapes(
    std::map<std::string, ShapeT *> &values, const char *marker)
{
    for (auto it = values.begin(); it != values.end();) {
	if (it->first.find(marker) == std::string::npos) {
	    ++it;
	    continue;
	}
	if (it->second)
	    it->second->unref();
	it = values.erase(it);
    }
}

void
BObolDatabaseSourceRealizationCache::eraseViewVariants(void)
{
    static const char marker[] = ":view-lod:";
    erase_realization_cache_variant_shapes(this->sharedWireGeometry, marker);
    erase_realization_cache_variant_values(this->sharedWireBounds, marker);
    erase_realization_cache_variant_shapes(this->sharedMeshVListGeometry,
	marker);
    erase_realization_cache_variant_shapes(this->sharedMeshGeometry, marker);
    erase_realization_cache_variant_values(this->sharedWireCadGeometry,
	marker);
    erase_realization_cache_variant_values(this->sharedMeshVListCadGeometry,
	marker);
    erase_realization_cache_variant_values(this->sharedMeshCadGeometry,
	marker);
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
    std::map<std::string, BObolCachedPartGeometry> &cache,
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
    std::map<std::string, BObolCachedPartGeometry> &cache,
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
    const std::map<std::string, BObolCachedPartGeometry> &cache,
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
