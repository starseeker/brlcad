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

#include <Obol/cad/SoCADAssembly.h>

#include <map>
#include <memory>
#include <string>

class SoBRLDatabaseSource;
class SoBRLVListShape;
class SoBRLMeshShape;
struct BObolCompactOccurrenceStream;

struct BObolCachedPartGeometry {
    BObolCachedPartGeometry(void) :
	geometryTransform(SbMatrix::identity()),
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
	const std::string &key, Obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL);
    std::shared_ptr<const Obol::PartGeometry> storeMeshVListCadGeometry(
	const std::string &key, Obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL);
    std::shared_ptr<const Obol::PartGeometry> storeMeshCadGeometry(
	const std::string &key, Obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL);
    std::shared_ptr<const Obol::PartGeometry> storeMeshCadGeometryReference(
	const std::string &key,
	const std::shared_ptr<const Obol::PartGeometry> &geometry,
	const SbMatrix &geometryTransform, const char *sourceType = NULL,
	const char *geometryKind = NULL, const SbBox3f *bounds = NULL,
	bool lodBacked = false,
	const BObolSourceMeshRequest *sourceMeshRequest = NULL);
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

#endif /* LIBBOBOL_DATABASE_SOURCE_REALIZATION_H */
