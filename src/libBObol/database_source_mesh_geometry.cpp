/*    D A T A B A S E _ S O U R C E _ M E S H _ G E O M E T R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLodRealization.h"
#include "BObol/BMeshLodCache.h"
#include "cad_normals_private.h"
#include "cad_publication_private.h"
#include "database_source_mesh_geometry_private.h"
#include "database_source_realization.h"
#include "performance_private.h"

#include "bg/trimesh.h"
#include "raytrace.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <unordered_set>
#include <utility>
#include <vector>

int
cad_wire_part_geometry_from_bot(const struct rt_bot_internal *bot,
	Obol::PartGeometryBuilder &geometry)
{
    if (!bot || !bot->vertices || !bot->faces || !bot->num_vertices ||
	!bot->num_faces ||
	bot->num_faces > std::numeric_limits<size_t>::max() / 6u ||
	bot->num_faces > std::numeric_limits<uint32_t>::max() / 3u)
	return 0;
    RT_BOT_CK_MAGIC(bot);

    Obol::WireRep wire;
    wire.bounds.makeEmpty();
    wire.segmentPoints.reserve(bot->num_faces * 6u);
    wire.segmentIds.reserve(bot->num_faces * 3u);
    uint32_t segmentId = 0;
    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    continue;
	const int edge[3][2] = {
	    {face[0], face[1]}, {face[1], face[2]}, {face[2], face[0]}
	};
	for (size_t j = 0; j < 3; j++) {
	    const fastf_t *a = &bot->vertices[edge[j][0] * 3];
	    const fastf_t *b = &bot->vertices[edge[j][1] * 3];
	    const SbVec3f start(static_cast<float>(a[X]),
		static_cast<float>(a[Y]), static_cast<float>(a[Z]));
	    const SbVec3f end(static_cast<float>(b[X]),
		static_cast<float>(b[Y]), static_cast<float>(b[Z]));
	    wire.segmentPoints.push_back(start);
	    wire.segmentPoints.push_back(end);
	    wire.segmentIds.push_back(segmentId++);
	    wire.bounds.extendBy(start);
	    wire.bounds.extendBy(end);
	}
    }
    if (wire.segmentPoints.empty())
	return 0;
    geometry.wire = std::move(wire);
    bobol_performance_counter_add(BOBOL_PERF_VLIST_POINTS,
	static_cast<uint64_t>(bot->num_faces) * 4u);
    return 1;
}

static bool
bot_authored_triangle_normals(std::vector<SbVec3f> &normals,
			     const struct rt_bot_internal *bot)
{
    normals.clear();
    if (!bot ||
	(bot->bot_flags &
	    (RT_BOT_HAS_SURFACE_NORMALS | RT_BOT_USE_NORMALS)) !=
	    (RT_BOT_HAS_SURFACE_NORMALS | RT_BOT_USE_NORMALS) ||
	!bot->normals || !bot->face_normals ||
	bot->num_face_normals < bot->num_faces)
	return false;

    normals.reserve(bot->num_faces * 3);
    for (size_t faceIndex = 0; faceIndex < bot->num_faces; faceIndex++) {
	for (size_t corner = 0; corner < 3; corner++) {
	    const size_t sourceCorner =
		(bot->orientation == RT_BOT_CW && corner > 0) ? 3 - corner : corner;
	    const int normalIndex =
		bot->face_normals[faceIndex * 3 + sourceCorner];
	    if (normalIndex < 0 ||
		static_cast<size_t>(normalIndex) >= bot->num_normals) {
		normals.clear();
		return false;
	    }
	    const fastf_t *normal = &bot->normals[static_cast<size_t>(normalIndex) * 3];
	    SbVec3f converted(static_cast<float>(normal[X]),
		static_cast<float>(normal[Y]), static_cast<float>(normal[Z]));
	    if (converted.length() <= 0.0f) {
		normals.clear();
		return false;
	    }
	    converted.normalize();
	    normals.push_back(converted);
	}
    }
    return true;
}

template <typename Index>
static void
bot_triangle_normals(std::vector<SbVec3f> &normals,
		     const struct rt_bot_internal *bot,
		     const std::vector<SbVec3f> &points,
		     const std::vector<Index> &triangles)
{
    if (!bot_authored_triangle_normals(normals, bot))
	generate_flat_triangle_normals(normals, points, triangles);
}

/* Keep the indexed type local to callers while sharing authored-normal
 * validation and winding correction. */
void
cad_bot_triangle_normals(std::vector<SbVec3f> &normals,
    const struct rt_bot_internal *bot,
    const std::vector<SbVec3f> &positions,
    const std::vector<int32_t> &indices)
{
    bot_triangle_normals(normals, bot, positions, indices);
}

void
cad_bot_triangle_normals(std::vector<SbVec3f> &normals,
    const struct rt_bot_internal *bot,
    const std::vector<SbVec3f> &positions,
    const std::vector<uint32_t> &indices)
{
    bot_triangle_normals(normals, bot, positions, indices);
}

int
cad_mesh_part_geometry_from_bot(const struct rt_bot_internal *bot,
	Obol::PartGeometryBuilder &geometry)
{
    if (!bot || !bot->vertices || !bot->faces ||
	bot->num_vertices == 0 || bot->num_faces == 0)
	return 0;
    RT_BOT_CK_MAGIC(bot);

    Obol::TriMesh mesh;
    mesh.bounds.makeEmpty();
    mesh.positions.reserve(bot->num_vertices);
    for (size_t i = 0; i < bot->num_vertices; i++) {
	const SbVec3f point(static_cast<float>(bot->vertices[i * 3]),
	    static_cast<float>(bot->vertices[i * 3 + 1]),
	    static_cast<float>(bot->vertices[i * 3 + 2]));
	mesh.positions.push_back(point);
	mesh.bounds.extendBy(point);
    }

    if (bot->num_faces > std::numeric_limits<size_t>::max() / 3)
	return 0;
    mesh.indices.reserve(bot->num_faces * 3);
    for (size_t i = 0; i < bot->num_faces; i++) {
	const int *face = &bot->faces[i * 3];
	if (face[0] < 0 || face[1] < 0 || face[2] < 0 ||
	    static_cast<size_t>(face[0]) >= bot->num_vertices ||
	    static_cast<size_t>(face[1]) >= bot->num_vertices ||
	    static_cast<size_t>(face[2]) >= bot->num_vertices)
	    return 0;
	mesh.indices.push_back(static_cast<uint32_t>(face[0]));
	if (bot->orientation == RT_BOT_CW) {
	    mesh.indices.push_back(static_cast<uint32_t>(face[2]));
	    mesh.indices.push_back(static_cast<uint32_t>(face[1]));
	} else {
	    mesh.indices.push_back(static_cast<uint32_t>(face[1]));
	    mesh.indices.push_back(static_cast<uint32_t>(face[2]));
	}
    }
    if (mesh.bounds.isEmpty() || mesh.indices.empty())
	return 0;
    bot_triangle_normals(mesh.normals, bot, mesh.positions, mesh.indices);
    if (!canonicalize_corner_normal_mesh(mesh, mesh.normals))
	return 0;
    geometry.shaded = std::move(mesh);
    geometry.shadedCullBackfaces =
	bot->mode == RT_BOT_SOLID &&
	bot->orientation != RT_BOT_UNORIENTED &&
	bot->num_vertices <= static_cast<size_t>(INT_MAX) &&
	bot->num_faces <= static_cast<size_t>(INT_MAX) &&
	bg_trimesh_solid2(static_cast<int>(bot->num_vertices),
	    static_cast<int>(bot->num_faces),
	    const_cast<fastf_t *>(bot->vertices),
	    const_cast<int *>(bot->faces), NULL) == 0;
    return 1;
}

int
cad_mesh_append_hidden_line_edges(Obol::PartGeometryBuilder &geometry)
{
    if (!geometry.shaded)
	return 0;
    const Obol::TriMesh &mesh = *geometry.shaded;
    Obol::WireRep wire;
    wire.bounds = mesh.bounds;
    std::unordered_set<uint64_t> edges;
    uint32_t segmentId = 0;
    auto addEdge = [&](uint32_t a, uint32_t b) {
	if (a >= mesh.positions.size() || b >= mesh.positions.size())
	    return;
	if (a > b)
	    std::swap(a, b);
	const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
	if (!edges.insert(key).second)
	    return;
	wire.segmentPoints.push_back(mesh.positions[a]);
	wire.segmentPoints.push_back(mesh.positions[b]);
	wire.segmentIds.push_back(segmentId++);
    };
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
	addEdge(mesh.indices[i], mesh.indices[i + 1]);
	addEdge(mesh.indices[i + 1], mesh.indices[i + 2]);
	addEdge(mesh.indices[i + 2], mesh.indices[i]);
    }
    if (wire.segmentPoints.empty())
	return 0;
    geometry.wire = std::move(wire);
    return 1;
}

static size_t
cad_geometry_saturating_add(size_t total, size_t count, size_t elementSize)
{
    if (!count || !elementSize || total == SIZE_MAX)
	return total;
    if (count > SIZE_MAX / elementSize)
	return SIZE_MAX;
    const size_t bytes = count * elementSize;
    return bytes > SIZE_MAX - total ? SIZE_MAX : total + bytes;
}

template <typename Geometry>
static size_t
cad_part_geometry_estimate_bytes(const Geometry &geometry)
{
    size_t bytes = sizeof(geometry);
    const auto addTriangleMesh = [&bytes](const Obol::TriMesh &mesh) {
	bytes = cad_geometry_saturating_add(bytes,
	    mesh.positions.capacity(), sizeof(SbVec3f));
	bytes = cad_geometry_saturating_add(bytes,
	    mesh.normals.capacity(), sizeof(SbVec3f));
	bytes = cad_geometry_saturating_add(bytes,
	    mesh.indices.capacity(), sizeof(uint32_t));
	bytes = cad_geometry_saturating_add(bytes,
	    mesh.progressiveCuts.capacity(),
	    sizeof(Obol::ProgressiveTriangleCut));
	bytes = cad_geometry_saturating_add(bytes,
	    mesh.progressiveClusters.capacity(),
	    sizeof(Obol::ProgressiveTriangleCluster));
	for (const Obol::ProgressiveTriangleCluster &cluster :
	     mesh.progressiveClusters)
	    bytes = cad_geometry_saturating_add(bytes,
		cluster.ranges.capacity(),
		sizeof(Obol::ProgressiveTriangleClusterRange));
    };

    if (geometry.points) {
	const Obol::PointRep &points = *geometry.points;
	bytes = cad_geometry_saturating_add(bytes,
	    points.positions.capacity(), sizeof(SbVec3f));
	bytes = cad_geometry_saturating_add(bytes,
	    points.pointIds.capacity(), sizeof(uint32_t));
	bytes = cad_geometry_saturating_add(bytes,
	    points.colorValid.capacity(), sizeof(uint8_t));
	bytes = cad_geometry_saturating_add(bytes,
	    points.colors.capacity(), sizeof(SbColor));
	bytes = cad_geometry_saturating_add(bytes,
	    points.scaleValid.capacity(), sizeof(uint8_t));
	bytes = cad_geometry_saturating_add(bytes,
	    points.scales.capacity(), sizeof(float));
	bytes = cad_geometry_saturating_add(bytes,
	    points.normalValid.capacity(), sizeof(uint8_t));
	bytes = cad_geometry_saturating_add(bytes,
	    points.normals.capacity(), sizeof(SbVec3f));
    }
    if (geometry.shaded)
	addTriangleMesh(*geometry.shaded);
    if (geometry.wire) {
	const Obol::WireRep &wire = *geometry.wire;
	bytes = cad_geometry_saturating_add(bytes,
	    wire.segmentPoints.capacity(), sizeof(SbVec3f));
	bytes = cad_geometry_saturating_add(bytes,
	    wire.segmentIds.capacity(), sizeof(uint32_t));
	bytes = cad_geometry_saturating_add(bytes,
	    wire.polylines.capacity(), sizeof(Obol::WirePolyline));
	for (const Obol::WirePolyline &polyline : wire.polylines)
	    bytes = cad_geometry_saturating_add(bytes,
		polyline.points.capacity(), sizeof(SbVec3f));
	bytes = cad_geometry_saturating_add(bytes,
	    wire.progressiveCuts.capacity(),
	    sizeof(Obol::ProgressiveWireCut));
	bytes = cad_geometry_saturating_add(bytes,
	    wire.progressiveClusters.capacity(),
	    sizeof(Obol::ProgressiveWireCluster));
	for (const Obol::ProgressiveWireCluster &cluster :
	     wire.progressiveClusters)
	    bytes = cad_geometry_saturating_add(bytes,
		cluster.ranges.capacity(),
		sizeof(Obol::ProgressiveWireClusterRange));
	/* A derived-wire source retains an independent admitted snapshot. */
	if (const Obol::TriMesh *triangleEdges = wire.triangleEdges())
	    addTriangleMesh(*triangleEdges);
    }
    return bytes;
}

size_t
bobol_database_part_geometry_estimate_bytes(
    const Obol::PartGeometry &geometry)
{
    return cad_part_geometry_estimate_bytes(geometry);
}

size_t
bobol_database_part_geometry_estimate_bytes(
    const Obol::PartGeometryBuilder &geometry)
{
    return cad_part_geometry_estimate_bytes(geometry);
}

std::shared_ptr<const Obol::PartGeometry>
bobol_database_bot_part_geometry(const struct rt_bot_internal *bot,
    int drawMode)
{
    Obol::PartGeometryBuilder geometry;
    int valid = 0;
    if (drawMode == BOBOL_LOD_DRAW_WIRE) {
	valid = cad_wire_part_geometry_from_bot(bot, geometry);
    } else if (drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE) {
	valid = cad_mesh_part_geometry_from_bot(bot, geometry);
	if (valid && drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	    valid = cad_mesh_append_hidden_line_edges(geometry);
    }
    if (!valid)
	return std::shared_ptr<const Obol::PartGeometry>();

    geometry.subpixelProxyEligible = true;
    return bobol_cad_build_geometry(
	std::move(geometry), "terminal BoT geometry");
}
