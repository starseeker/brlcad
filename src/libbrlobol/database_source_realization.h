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

#include <obol/cad/SoCADAssembly.h>

#include <map>
#include <string>

class SoBRLDatabaseSource;
class SoBRLVListShape;
class SoBRLMeshShape;

struct BRLObolDatabaseSourceRealizationCache {
    BRLObolDatabaseSourceRealizationCache(void);
    ~BRLObolDatabaseSourceRealizationCache(void);

    void clear(void);
    void storeWireGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeWireBounds(const std::string &key, const SbBox3f &bounds);
    void storeMeshVListGeometry(const std::string &key, SoBRLVListShape *shape);
    void storeMeshGeometry(const std::string &key, SoBRLMeshShape *shape);
    void storeWireCadGeometry(const std::string &key,
	const obol::PartGeometry &geometry);

    std::map<std::string, SoBRLVListShape *> sharedWireGeometry;
    std::map<std::string, SbBox3f> sharedWireBounds;
    std::map<std::string, SoBRLVListShape *> sharedMeshVListGeometry;
    std::map<std::string, SoBRLMeshShape *> sharedMeshGeometry;
    std::map<std::string, obol::PartGeometry> sharedWireCadGeometry;
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
