/*        D A T A B A S E _ S O U R C E _ R E A L I Z A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DATABASE_SOURCE_REALIZATION_H
#define LIBBOBOL_DATABASE_SOURCE_REALIZATION_H

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbMatrix.h>

#include "BObol/BSourceMeshRequest.h"

#include "bg/defines.h"

#include <Obol/cad/SoCADAssembly.h>

#include <map>
#include <exception>
#include <memory>
#include <string>
#include <vector>

class SbString;
class SoNode;
class SoGroup;
class SoBRLDatabaseSource;
class SoBRLVListShape;
class SoBRLMeshShape;
struct BObolCompactOccurrenceStream;
struct BObolCompactOccurrence;
struct rt_bot_internal;
struct db_i;
struct BObolSourceRealizationCoordinatorPrivate;

/* Private stream-to-queue cancellation edge. The weak endpoint cannot keep
 * the process-wide pool alive; cleanup runs outside both owner mutexes. */
void bobol_source_realization_cancel_queued(
    const std::weak_ptr<BObolSourceRealizationCoordinatorPrivate> &coordinator);

/* Streamed v5 BoTs use serialized coverage and mesh-service admission.  The
 * source worker must not enter the direct primitive-import fallback. */
bool bobol_database_source_uses_serialized_bot_coverage(
    const SoBRLDatabaseSource *source, struct db_i *database);

/* Whole-target coverage bypasses the leaf subpixel classifier and can be
 * delivered independently of storage for the full occurrence population. */
BObolCompactOccurrence bobol_database_source_coverage_overview(
    const SoBRLDatabaseSource *source, const char *treeName,
    const SbBox3f &bounds, uint32_t revision);

struct BObolCachedPartGeometry {
    BObolCachedPartGeometry(void) :
	geometryTransform(SbMatrix::identity()),
	viewDependentCsgGeometry(false),
	lodBacked(false),
	sourceMeshRequestValid(false)
    {
	bounds.makeEmpty();
    }

    std::shared_ptr<const Obol::PartGeometry> geometry;
    std::string sourceType;
    std::string geometryKind;
    SbBox3f bounds;
    /* Maps the shared geometry's local coordinates into this cache key's
     * object-local coordinates.  Identity for ordinary cache entries. */
    SbMatrix geometryTransform;
    bool viewDependentCsgGeometry;
    bool lodBacked;
    bool sourceMeshRequestValid;
    BObolSourceMeshRequest sourceMeshRequest;
};

template <typename Value>
using BObolRealizationCacheMap = std::map<std::string, Value, std::less<>>;

/* Cache retention must not rewrite an already published geometry's identity.
 * Producers initialize identity before storing; seeding only shares ownership. */
template <typename Shape>
void bobol_cache_geometry_reference(BObolRealizationCacheMap<Shape *> &cache,
    const std::string &key, Shape *shape)
{
    if (!shape) return;
    const auto found = cache.find(key);
    Shape *previous = found == cache.end() ? nullptr : found->second;
    if (previous == shape) return;
    shape->ref();
    try { cache.insert_or_assign(key, shape); }
    catch (...) { shape->unref(); throw; }
    if (previous) previous->unref();
}

struct BObolDatabaseSourceRealizationCache {
    BObolDatabaseSourceRealizationCache(void);
    ~BObolDatabaseSourceRealizationCache(void);
    BObolDatabaseSourceRealizationCache(
	const BObolDatabaseSourceRealizationCache &) = delete;
    BObolDatabaseSourceRealizationCache &operator=(
	const BObolDatabaseSourceRealizationCache &) = delete;

    void clear(void) noexcept;
    void eraseObject(const std::string &name);
    void renameObject(const std::string &oldName,
	const std::string &newName);
    class ObjectRename {
    public:
	~ObjectRename();
	void commit() noexcept;
    private:
	friend struct BObolDatabaseSourceRealizationCache;
	struct Impl;
	explicit ObjectRename(std::unique_ptr<Impl> prepared);
	std::unique_ptr<Impl> d;
    };
    std::unique_ptr<ObjectRename> prepareObjectRename(
	const std::string &oldName, const std::string &newName);
    void eraseViewVariants(void) noexcept;
    void swap(BObolDatabaseSourceRealizationCache &other) noexcept;
    void storeWireGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeWireBounds(const std::string &key, const SbBox3f &bounds);
    void storeMeshVListGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeMeshGeometry(const std::string &key, SoBRLMeshShape *shape);
    std::shared_ptr<const Obol::PartGeometry> storeWireCadGeometry(
	const std::string &key, Obol::PartGeometryBuilder &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL,
	bool viewDependentCsgGeometry = false);
    std::shared_ptr<const Obol::PartGeometry> storeMeshVListCadGeometry(
	const std::string &key, Obol::PartGeometryBuilder &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL,
	bool viewDependentCsgGeometry = false);
    std::shared_ptr<const Obol::PartGeometry> storeMeshCadGeometry(
	const std::string &key, Obol::PartGeometryBuilder &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL,
	bool viewDependentCsgGeometry = false);
    std::shared_ptr<const Obol::PartGeometry> storeMeshCadGeometryReference(
	const std::string &key,
	const std::shared_ptr<const Obol::PartGeometry> &geometry,
	const SbMatrix &geometryTransform, const char *sourceType = NULL,
	const char *geometryKind = NULL, const SbBox3f *bounds = NULL,
	bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL,
	bool viewDependentCsgGeometry = false);
    const BObolCachedPartGeometry *findWireCadGeometry(
	const std::string &key) const;
    const BObolCachedPartGeometry *findMeshVListCadGeometry(
	const std::string &key) const;
    const BObolCachedPartGeometry *findMeshCadGeometry(
	const std::string &key) const;

    BObolRealizationCacheMap<SoBRLVListShape *> sharedWireGeometry;
    BObolRealizationCacheMap<SbBox3f> sharedWireBounds;
    BObolRealizationCacheMap<SoBRLVListShape *> sharedMeshVListGeometry;
    BObolRealizationCacheMap<SoBRLMeshShape *> sharedMeshGeometry;
    BObolRealizationCacheMap<BObolCachedPartGeometry>
	sharedWireCadGeometry;
    BObolRealizationCacheMap<BObolCachedPartGeometry>
	sharedMeshVListCadGeometry;
    BObolRealizationCacheMap<BObolCachedPartGeometry>
	sharedMeshCadGeometry;
};

std::string bobol_realization_cache_object_name(const std::string &path);

/* Exclusive private construction: source has no live observers and is owned
 * by the producer until terminal handoff. Public realization APIs instead
 * prepare and publish a complete replacement for an observable source. */
SbBool bobol_database_source_construct_realization(SoBRLDatabaseSource *source,
	SbBool mesh, BObolCompactOccurrenceStream *stream);

/* Borrowed for one synchronous owner-thread publication. Successful
 * realization supplies the exact prepared child order and the retained nodes
 * whose metadata changes. Preparation can fail; commit runs after all source
 * writes and before any observer notification. Failed realization supplies
 * empty child and changed-node lists. */
struct BObolSourceRealizationEffects {
    virtual ~BObolSourceRealizationEffects() = default;
    virtual void prepare(const SoBRLDatabaseSource &source, bool realized,
	const SbString &diagnostic, const std::vector<SoNode *> &children,
	const std::vector<SoNode *> &changedNodes) = 0;
    virtual void commit(bool changed) noexcept = 0;
    virtual void notify() = 0;
    void notify(std::exception_ptr &failure) noexcept
    {
	try { this->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
};

/* Source child writers stage every affected group before changing the live
 * graph. The scene uses the complete set to prepare indexes, repository
 * membership and peer revisions, then commits those effects before observers. */
struct BObolSourceChildEffects {
    virtual ~BObolSourceChildEffects() = default;
    virtual void stageChildOrder(SoGroup &parent,
	const std::vector<SoNode *> &children) = 0;
    virtual void stageFrameEffect(SoNode &node) = 0;
    virtual void prepare() = 0;
    virtual void commit() noexcept = 0;
};

/* A resident-occurrence adoption prepares its enclosing scene retirement
 * before changing the target registry. The source commits the registry first,
 * then this effect retires donor edges and advances the scene before either
 * source or hierarchy observers run. */
struct BObolSourceAdoptionEffects {
    virtual ~BObolSourceAdoptionEffects() = default;
    virtual void prepare() = 0;
    virtual void commit() noexcept = 0;
    virtual void notify() = 0;
    void notify(std::exception_ptr &failure) noexcept
    {
	try { this->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
};

SbBool bobol_database_source_realize_prototype(
	SoBRLDatabaseSource *source,
	BObolSourceRealizationEffects *effects = nullptr);
SbBool bobol_database_source_realize_wireframe_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolSourceRealizationEffects *effects = nullptr);
SbBool bobol_database_source_realize_mesh_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolSourceRealizationEffects *effects = nullptr);
int bobol_database_source_realize_wireframe_compact_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolCompactOccurrenceStream *stream = NULL,
	BObolSourceRealizationEffects *effects = nullptr);
int bobol_database_source_realize_mesh_compact_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolCompactOccurrenceStream *stream = NULL,
	BObolSourceRealizationEffects *effects = nullptr);
void bobol_database_source_seed_realization_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache);

/* Build the immutable renderer representation of one terminal BoT on a
 * worker.  Keeping this conversion beside ordinary database realization
 * ensures winding, normals, hidden-line edges, and culling certification use
 * one implementation. */
std::shared_ptr<const Obol::PartGeometry>
bobol_database_bot_part_geometry(const struct rt_bot_internal *bot,
	int drawMode);
size_t bobol_database_part_geometry_estimate_bytes(
	const Obol::PartGeometry &geometry);
size_t bobol_database_part_geometry_estimate_bytes(
	const Obol::PartGeometryBuilder &geometry);

/* Generate one detached BREP triangle representation.  The caller supplies
 * its deterministic band identity; the returned owner releases all
 * tessellation arrays after the PoP cache has consumed them. */
std::shared_ptr<BObolStagedSourceMesh>
bobol_database_brep_staged_mesh_variant(
	struct db_i *dbip, const char *name,
	const struct bg_tess_tol *ttol, uint64_t contentKey,
	uint32_t sourceRevision, BObolSourceMeshRequest &request);

#endif /* LIBBOBOL_DATABASE_SOURCE_REALIZATION_H */
