/* D A T A B A S E _ S O U R C E _ M E S H _ G E O M E T R Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_DATABASE_SOURCE_MESH_GEOMETRY_PRIVATE_H
#define LIBBOBOL_DATABASE_SOURCE_MESH_GEOMETRY_PRIVATE_H

#include <Inventor/SbVec3f.h>
#include <Obol/cad/CadGeometry.h>

#include <cstdint>
#include <vector>

struct rt_bot_internal;

int cad_wire_part_geometry_from_bot(const struct rt_bot_internal *bot,
    Obol::PartGeometryBuilder &geometry);
int cad_mesh_part_geometry_from_bot(const struct rt_bot_internal *bot,
    Obol::PartGeometryBuilder &geometry);
int cad_mesh_append_hidden_line_edges(Obol::PartGeometryBuilder &geometry);

void cad_bot_triangle_normals(std::vector<SbVec3f> &normals,
    const struct rt_bot_internal *bot,
    const std::vector<SbVec3f> &positions,
    const std::vector<int32_t> &indices);
void cad_bot_triangle_normals(std::vector<SbVec3f> &normals,
    const struct rt_bot_internal *bot,
    const std::vector<SbVec3f> &positions,
    const std::vector<uint32_t> &indices);

#endif /* LIBBOBOL_DATABASE_SOURCE_MESH_GEOMETRY_PRIVATE_H */
