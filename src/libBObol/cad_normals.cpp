/*                    C A D _ N O R M A L S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "cad_normals_private.h"

#include <cstring>
#include <map>
#include <utility>

void
sanitize_triangle_normals(std::vector<SbVec3f> &normals,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &triangles)
{
    if (normals.empty()) {
	generate_flat_triangle_normals(normals, points, triangles);
	return;
    }
    if (normals.size() != triangles.size())
	return;

    for (size_t i = 0; i + 2 < triangles.size(); i += 3) {
	const int32_t ia = triangles[i];
	const int32_t ib = triangles[i + 1];
	const int32_t ic = triangles[i + 2];
	if (ia < 0 || ib < 0 || ic < 0 ||
	    static_cast<size_t>(ia) >= points.size() ||
	    static_cast<size_t>(ib) >= points.size() ||
	    static_cast<size_t>(ic) >= points.size())
	    continue;

	SbVec3f faceNormal =
	    (points[static_cast<size_t>(ib)] - points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] - points[static_cast<size_t>(ia)]);
	if (faceNormal.length() > 0.0f)
	    faceNormal.normalize();
	else
	    continue;

	float normalDot = 0.0f;
	for (size_t j = 0; j < 3; j++) {
	    SbVec3f &n = normals[i + j];
	    if (n.length() > 0.0f) {
		n.normalize();
		normalDot += n.dot(faceNormal);
	}
}

	if (normalDot < 0.0f) {
	    for (size_t j = 0; j < 3; j++)
		normals[i + j].negate();
	}

	for (size_t j = 0; j < 3; j++) {
	    if (normals[i + j].length() <= 0.0f)
		normals[i + j] = faceNormal;
	}
    }
}

static uint32_t
corner_normal_float_bits(float value)
{
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

struct CornerNormalVertexKey {
    uint32_t position;
    uint32_t normalX;
    uint32_t normalY;
    uint32_t normalZ;

    bool operator<(const CornerNormalVertexKey &other) const
    {
	if (position != other.position)
	    return position < other.position;
	if (normalX != other.normalX)
	    return normalX < other.normalX;
	if (normalY != other.normalY)
	    return normalY < other.normalY;
	return normalZ < other.normalZ;
    }
};

/* Obol's indexed mesh format binds one normal to each vertex.  BRL-CAD mesh
 * producers keep normals per triangle corner, so split only corners whose
 * normal differs instead of losing smooth shading or duplicating every vertex. */
int
canonicalize_corner_normal_mesh(Obol::TriMesh &mesh,
	const std::vector<SbVec3f> &cornerNormals)
{
    if (cornerNormals.empty()) {
	mesh.normals.clear();
	return 1;
    }
    if (cornerNormals.size() != mesh.indices.size())
	return 0;

    std::map<CornerNormalVertexKey, uint32_t> vertexByCorner;
    std::vector<SbVec3f> positions;
    std::vector<SbVec3f> normals;
    std::vector<uint32_t> indices;
    positions.reserve(mesh.indices.size());
    normals.reserve(mesh.indices.size());
    indices.reserve(mesh.indices.size());
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
	const uint32_t sourceIndex = mesh.indices[i];
	if (sourceIndex >= mesh.positions.size())
	    return 0;
	SbVec3f normal = cornerNormals[i];
	if (normal.sqrLength() > 0.0f)
	    normal.normalize();
	else
	    normal.setValue(0.0f, 0.0f, 1.0f);
	const CornerNormalVertexKey key = {
	    sourceIndex,
	    corner_normal_float_bits(normal[0]),
	    corner_normal_float_bits(normal[1]),
	    corner_normal_float_bits(normal[2])
	};
	auto found = vertexByCorner.find(key);
	if (found == vertexByCorner.end()) {
	    const uint32_t index = static_cast<uint32_t>(positions.size());
	    positions.push_back(mesh.positions[sourceIndex]);
	    normals.push_back(normal);
	    found = vertexByCorner.emplace(key, index).first;
	}
	indices.push_back(found->second);
    }

    mesh.positions = std::move(positions);
    mesh.normals = std::move(normals);
    mesh.indices = std::move(indices);
    return 1;
}
