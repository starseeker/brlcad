/* C A D _ N O R M A L S _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_CAD_NORMALS_PRIVATE_H
#define LIBBOBOL_CAD_NORMALS_PRIVATE_H

#include <Inventor/SbVec3f.h>
#include <Obol/cad/CadGeometry.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

void sanitize_triangle_normals(std::vector<SbVec3f> &normals,
    const std::vector<SbVec3f> &points,
    const std::vector<int32_t> &triangles);

/* Convert per-corner normals into Obol's one-normal-per-indexed-vertex
 * representation, splitting only position/normal pairs which differ. */
int canonicalize_corner_normal_mesh(Obol::TriMesh &mesh,
    const std::vector<SbVec3f> &cornerNormals);

/* Missing normals carry no author intent to smooth.  Preserve BRL-CAD's
 * default shaded-BoT semantics by using the geometric normal of each triangle
 * at all three corners.  The separate per-view smooth policy may explicitly
 * synthesize presentation normals later. */
template <typename Index>
static void
generate_flat_triangle_normals(std::vector<SbVec3f> &normals,
			       const std::vector<SbVec3f> &points,
			       const std::vector<Index> &triangles)
{
    normals.clear();
    if (points.empty() || triangles.empty() || triangles.size() % 3)
	return;

    const size_t faceCount = triangles.size() / 3;
    normals.reserve(triangles.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; faceIndex++) {
	const size_t indexBase = faceIndex * 3;
	const uint64_t ia = static_cast<uint64_t>(triangles[indexBase]);
	const uint64_t ib = static_cast<uint64_t>(triangles[indexBase + 1]);
	const uint64_t ic = static_cast<uint64_t>(triangles[indexBase + 2]);
	if (ia >= points.size() || ib >= points.size() || ic >= points.size()) {
	    normals.clear();
	    return;
	}

	SbVec3f normal =
	    (points[static_cast<size_t>(ib)] - points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] - points[static_cast<size_t>(ia)]);
	if (normal.length() > 0.0f)
	    normal.normalize();
	else
	    normal = SbVec3f(0.0f, 0.0f, 1.0f);
	normals.push_back(normal);
	normals.push_back(normal);
	normals.push_back(normal);
    }
}

template <typename Index>
static void
generate_smooth_triangle_normals(std::vector<SbVec3f> &normals,
				 const std::vector<SbVec3f> &points,
				 const std::vector<Index> &triangles,
				 float creaseAngleDegrees)
{
    normals.clear();
    if (points.empty() || triangles.empty() || triangles.size() % 3)
	return;

    const size_t faceCount = triangles.size() / 3;
    std::vector<SbVec3f> faceNormals(faceCount);
    std::vector<std::vector<size_t>> vertexFaces(points.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const size_t base = faceIndex * 3;
	const uint64_t ia = static_cast<uint64_t>(triangles[base]);
	const uint64_t ib = static_cast<uint64_t>(triangles[base + 1]);
	const uint64_t ic = static_cast<uint64_t>(triangles[base + 2]);
	if (ia >= points.size() || ib >= points.size() || ic >= points.size()) {
	    normals.clear();
	    return;
	}
	SbVec3f normal =
	    (points[static_cast<size_t>(ib)] -
	     points[static_cast<size_t>(ia)]).cross(
		points[static_cast<size_t>(ic)] -
		points[static_cast<size_t>(ia)]);
	if (normal.sqrLength() > 0.0f)
	    normal.normalize();
	else
	    normal.setValue(0.0f, 0.0f, 1.0f);
	faceNormals[faceIndex] = normal;
	vertexFaces[static_cast<size_t>(ia)].push_back(faceIndex);
	vertexFaces[static_cast<size_t>(ib)].push_back(faceIndex);
	vertexFaces[static_cast<size_t>(ic)].push_back(faceIndex);
    }

    const float clampedAngle =
	std::max(0.0f, std::min(180.0f, creaseAngleDegrees));
    const float creaseCosine =
	static_cast<float>(std::cos(clampedAngle * M_PI / 180.0));
    normals.reserve(triangles.size());
    for (size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex) {
	const SbVec3f &faceNormal = faceNormals[faceIndex];
	for (size_t corner = 0; corner < 3; ++corner) {
	    const size_t vertexIndex = static_cast<size_t>(
		triangles[faceIndex * 3 + corner]);
	    SbVec3f smoothNormal(0.0f, 0.0f, 0.0f);
	    for (const size_t adjacentFace : vertexFaces[vertexIndex]) {
		const SbVec3f &candidate = faceNormals[adjacentFace];
		if (candidate.dot(faceNormal) >= creaseCosine)
		    smoothNormal += candidate;
	    }
	    if (smoothNormal.sqrLength() > 0.0f)
		smoothNormal.normalize();
	    else
		smoothNormal = faceNormal;
	    normals.push_back(smoothNormal);
	}
    }
}

#endif /* LIBBOBOL_CAD_NORMALS_PRIVATE_H */
