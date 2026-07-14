/*        D A T A B A S E _ S O U R C E _ R E A L I Z A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBRLOBOL_DATABASE_SOURCE_REALIZATION_H
#define LIBBRLOBOL_DATABASE_SOURCE_REALIZATION_H

#include <Inventor/SbBasic.h>
#include <Inventor/SbBox.h>

#include "brlobol/source_mesh_request.h"

#include <obol/cad/SoCADAssembly.h>

#include <map>
#include <memory>
#include <string>

class SoBRLDatabaseSource;
class SoBRLVListShape;
class SoBRLMeshShape;

struct BRLObolCachedPartGeometry {
    BRLObolCachedPartGeometry(void) :
	lodBacked(false),
	sourceMeshRequestValid(false)
    {
	bounds.makeEmpty();
    }

    std::shared_ptr<const obol::PartGeometry> geometry;
    std::string sourceType;
    std::string geometryKind;
    SbBox3f bounds;
    bool lodBacked;
    bool sourceMeshRequestValid;
    BRLObolSourceMeshRequest sourceMeshRequest;
};

struct BRLObolDatabaseSourceRealizationCache {
    BRLObolDatabaseSourceRealizationCache(void);
    ~BRLObolDatabaseSourceRealizationCache(void);

    void clear(void);
    void eraseObject(const std::string &name);
    void renameObject(const std::string &oldName,
	const std::string &newName);
    void eraseViewVariants(void);
    void storeWireGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeWireBounds(const std::string &key, const SbBox3f &bounds);
    void storeMeshVListGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeMeshGeometry(const std::string &key, SoBRLMeshShape *shape);
    std::shared_ptr<const obol::PartGeometry> storeWireCadGeometry(
	const std::string &key, obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BRLObolSourceMeshRequest *sourceMeshRequest = NULL);
    std::shared_ptr<const obol::PartGeometry> storeMeshVListCadGeometry(
	const std::string &key, obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BRLObolSourceMeshRequest *sourceMeshRequest = NULL);
    std::shared_ptr<const obol::PartGeometry> storeMeshCadGeometry(
	const std::string &key, obol::PartGeometry &&geometry,
	const char *sourceType = NULL, const char *geometryKind = NULL,
	const SbBox3f *bounds = NULL, bool lodBacked = false,
	const BRLObolSourceMeshRequest *sourceMeshRequest = NULL);
    const BRLObolCachedPartGeometry *findWireCadGeometry(
	const std::string &key) const;
    const BRLObolCachedPartGeometry *findMeshVListCadGeometry(
	const std::string &key) const;
    const BRLObolCachedPartGeometry *findMeshCadGeometry(
	const std::string &key) const;

    std::map<std::string, SoBRLVListShape *> sharedWireGeometry;
    std::map<std::string, SbBox3f> sharedWireBounds;
    std::map<std::string, SoBRLVListShape *> sharedMeshVListGeometry;
    std::map<std::string, SoBRLMeshShape *> sharedMeshGeometry;
    std::map<std::string, BRLObolCachedPartGeometry>
	sharedWireCadGeometry;
    std::map<std::string, BRLObolCachedPartGeometry>
	sharedMeshVListCadGeometry;
    std::map<std::string, BRLObolCachedPartGeometry>
	sharedMeshCadGeometry;
    bool preserveCompactSourceOnFailure;
};

SbBool brlobol_database_source_realize_wireframe_with_cache(
	SoBRLDatabaseSource *source,
	BRLObolDatabaseSourceRealizationCache *cache);
SbBool brlobol_database_source_realize_mesh_with_cache(
	SoBRLDatabaseSource *source,
	BRLObolDatabaseSourceRealizationCache *cache);
int brlobol_database_source_realize_wireframe_compact_with_cache(
	SoBRLDatabaseSource *source,
	BRLObolDatabaseSourceRealizationCache *cache);
int brlobol_database_source_realize_mesh_compact_with_cache(
	SoBRLDatabaseSource *source,
	BRLObolDatabaseSourceRealizationCache *cache);
void brlobol_database_source_seed_realization_cache(
	SoBRLDatabaseSource *source,
	BRLObolDatabaseSourceRealizationCache *cache);

#endif /* LIBBRLOBOL_DATABASE_SOURCE_REALIZATION_H */
