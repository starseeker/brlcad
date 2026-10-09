/*       C A D _ P U B L I C A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_CAD_PUBLICATION_PRIVATE_H
#define LIBBOBOL_CAD_PUBLICATION_PRIVATE_H

#include <Obol/cad/CadGeometryValidation.h>
#include <Obol/cad/CadSceneMutation.h>
#include <Obol/cad/CadSceneRecords.h>

#include <functional>
#include <memory>
#include <vector>

namespace Obol {
class PartGeometry;
struct PartGeometryBuilder;
}

class SoCADAssembly;
struct BObolLodCounts;

bool bobol_cad_admit_geometry(
    const std::shared_ptr<const Obol::PartGeometry> &geometry,
    Obol::ValidatedPartGeometry &admitted, const char *operation);

std::shared_ptr<const Obol::PartGeometry> bobol_cad_build_geometry(
    Obol::PartGeometryBuilder geometry, const char *operation);

std::shared_ptr<const Obol::PartGeometry>
bobol_cad_build_geometry_with_optional_proxy(
    Obol::PartGeometryBuilder geometry, const char *operation);

/** Describe the primitive/vertex work submitted by one immutable part. */
BObolLodCounts bobol_cad_geometry_counts(
    const Obol::PartGeometry &geometry);

bool bobol_cad_validate_shared_parts(
    const std::vector<Obol::PartUpdate> &updates,
    const char *operation);

bool bobol_cad_validate_instances(
    const std::vector<Obol::InstanceUpdate> &updates,
    const char *operation);

bool bobol_cad_validate_instance(const Obol::InstanceUpdate &update,
    size_t updateIndex, const char *operation);

bool bobol_cad_validate_styles(
    const std::vector<Obol::InstanceStyleUpdate> &updates,
    const char *operation);

bool bobol_cad_publish_styles(SoCADAssembly *assembly,
    const std::vector<Obol::InstanceStyleUpdate> &updates,
    const char *operation);

bool bobol_cad_publish_mutation(SoCADAssembly *assembly,
    const Obol::CadSceneMutation &mutation, const char *operation);

bool bobol_cad_validate_mutation(const SoCADAssembly *assembly,
    const Obol::CadSceneMutation &mutation, const char *operation);

bool bobol_cad_replace_scene(SoCADAssembly *assembly,
    const std::vector<Obol::PartUpdate> &parts,
    const std::vector<Obol::InstanceUpdate> &instances,
    const char *operation);

bool bobol_cad_replace_scene(SoCADAssembly *assembly,
    const std::vector<Obol::PartUpdate> &parts, size_t instanceCount,
    const std::function<Obol::InstanceUpdate(size_t)> &instanceAt,
    const char *operation);

#endif /* LIBBOBOL_CAD_PUBLICATION_PRIVATE_H */
