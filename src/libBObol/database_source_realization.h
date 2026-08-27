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
#include <memory>
#include <string>

class SoBRLDatabaseSource;
class SoBRLVListShape;
class SoBRLMeshShape;
struct BObolCompactOccurrenceStream;
struct rt_bot_internal;

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

struct BObolDatabaseSourceRealizationCache {
    BObolDatabaseSourceRealizationCache(void);
    ~BObolDatabaseSourceRealizationCache(void);

    void clear(void);
    void eraseObject(const std::string &name);
    void renameObject(const std::string &oldName,
	const std::string &newName);
    void eraseViewVariants(void);
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

    std::map<std::string, SoBRLVListShape *> sharedWireGeometry;
    std::map<std::string, SbBox3f> sharedWireBounds;
    std::map<std::string, SoBRLVListShape *> sharedMeshVListGeometry;
    std::map<std::string, SoBRLMeshShape *> sharedMeshGeometry;
    std::map<std::string, BObolCachedPartGeometry>
	sharedWireCadGeometry;
    std::map<std::string, BObolCachedPartGeometry>
	sharedMeshVListCadGeometry;
    std::map<std::string, BObolCachedPartGeometry>
	sharedMeshCadGeometry;
    bool preserveCompactSourceOnFailure;
};

SbBool bobol_database_source_realize_wireframe_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache);
SbBool bobol_database_source_realize_mesh_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache);
int bobol_database_source_realize_wireframe_compact_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolCompactOccurrenceStream *stream = NULL);
int bobol_database_source_realize_mesh_compact_with_cache(
	SoBRLDatabaseSource *source,
	BObolDatabaseSourceRealizationCache *cache,
	BObolCompactOccurrenceStream *stream = NULL);
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
