/*                 T R I M E S H _ D I F F . C P P
 * BRL-CAD
 *
 * Copyright (c) 2025-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file trimesh_diff.cpp
 *
 * Given two arrays of faces and points, determine if they define
 * the same mesh with tolerance.
 */

#include "common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

#include "bu/hash.h"
#include "bg/pca.h"
#include "bg/trimesh.h"

namespace {

struct tm_diff_vert {
    fastf_t x;
    fastf_t y;
    fastf_t z;
    size_t orig_ind;
};

typedef std::array<size_t, 3> tm_triangle;

static bool
vertex_less(const tm_diff_vert &a, const tm_diff_vert &b)
{
    const std::less<fastf_t> less;
    if (less(a.x, b.x))
	return true;
    if (less(b.x, a.x))
	return false;
    if (less(a.y, b.y))
	return true;
    if (less(b.y, a.y))
	return false;
    if (less(a.z, b.z))
	return true;
    if (less(b.z, a.z))
	return false;
    return a.orig_ind < b.orig_ind;
}

static bool
vertex_equal(const tm_diff_vert &a, const tm_diff_vert &b, fastf_t tol)
{
    return NEAR_EQUAL(a.x, b.x, tol) &&
	NEAR_EQUAL(a.y, b.y, tol) &&
	NEAR_EQUAL(a.z, b.z, tol);
}

static bool
valid_points(const point_t *points, size_t pointCount)
{
    if (!pointCount)
	return true;
    if (!points)
	return false;
    for (size_t i = 0; i < pointCount; i++) {
	if (!std::isfinite(points[i][X]) || !std::isfinite(points[i][Y]) ||
	    !std::isfinite(points[i][Z]))
	    return false;
    }
    return true;
}

static bool
valid_faces(const int *faces, size_t faceCount, size_t pointCount)
{
    if (!faceCount)
	return true;
    if (!faces || !pointCount)
	return false;
    if (faceCount > std::numeric_limits<size_t>::max() / 3u)
	return false;
    for (size_t i = 0; i < faceCount * 3u; i++) {
	if (faces[i] < 0 || static_cast<size_t>(faces[i]) >= pointCount)
	    return false;
    }
    return true;
}

static void
sorted_vertices(std::vector<tm_diff_vert> &vertices, const point_t *points,
	size_t pointCount)
{
    vertices.clear();
    vertices.reserve(pointCount);
    for (size_t i = 0; i < pointCount; i++) {
	tm_diff_vert vertex;
	vertex.x = points[i][X];
	vertex.y = points[i][Y];
	vertex.z = points[i][Z];
	vertex.orig_ind = i;
	vertices.push_back(vertex);
    }
    std::sort(vertices.begin(), vertices.end(), vertex_less);
}

static void
rank_by_original_index(std::vector<size_t> &ranks,
	const std::vector<tm_diff_vert> &vertices)
{
    ranks.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++)
	ranks[vertices[i].orig_ind] = i;
}

static bool
points_match_by_index(const point_t *first, const point_t *second,
	size_t pointCount, fastf_t tolerance)
{
    for (size_t i = 0; i < pointCount; i++) {
	tm_diff_vert a;
	tm_diff_vert b;
	a.x = first[i][X];
	a.y = first[i][Y];
	a.z = first[i][Z];
	b.x = second[i][X];
	b.y = second[i][Y];
	b.z = second[i][Z];
	if (!vertex_equal(a, b, tolerance))
	    return false;
    }
    return true;
}

static tm_triangle
canonical_triangle(size_t a, size_t b, size_t c)
{
    const tm_triangle candidates[3] = {{a, b, c}, {b, c, a}, {c, a, b}};
    return std::min(std::min(candidates[0], candidates[1]), candidates[2]);
}

static bool
canonical_faces(std::vector<tm_triangle> &triangles, const int *faces,
	size_t faceCount, const std::vector<size_t> &ranks)
{
    triangles.clear();
    triangles.reserve(faceCount);
    for (size_t i = 0; i < faceCount; i++) {
	const size_t first = static_cast<size_t>(faces[i * 3u]);
	const size_t second = static_cast<size_t>(faces[i * 3u + 1u]);
	const size_t third = static_cast<size_t>(faces[i * 3u + 2u]);
	if (first >= ranks.size() || second >= ranks.size() ||
	    third >= ranks.size())
	    return false;
	triangles.push_back(canonical_triangle(ranks[first], ranks[second],
	    ranks[third]));
    }
    std::sort(triangles.begin(), triangles.end());
    return true;
}

static bool
quantize_coordinate(fastf_t value, fastf_t tolerance, int64_t &quantized)
{
    const long double scaled = static_cast<long double>(value) /
	static_cast<long double>(tolerance);
    if (!std::isfinite(scaled) ||
	scaled > static_cast<long double>(std::numeric_limits<int64_t>::max()) ||
	scaled < static_cast<long double>(std::numeric_limits<int64_t>::min()))
	return false;
    quantized = static_cast<int64_t>(std::round(scaled));
    return true;
}

struct tm_hash_vert {
    int64_t x;
    int64_t y;
    int64_t z;
    size_t orig_ind;
};

struct tm_quantized_point {
    int64_t x;
    int64_t y;
    int64_t z;
};

static bool
hash_vertex_less(const tm_hash_vert &a, const tm_hash_vert &b)
{
    if (a.x != b.x)
	return a.x < b.x;
    if (a.y != b.y)
	return a.y < b.y;
    if (a.z != b.z)
	return a.z < b.z;
    return a.orig_ind < b.orig_ind;
}

static bool
hash_vertices(std::vector<tm_hash_vert> &vertices, const point_t *points,
	size_t pointCount, fastf_t tolerance)
{
    vertices.clear();
    vertices.reserve(pointCount);
    for (size_t i = 0; i < pointCount; i++) {
	tm_hash_vert vertex;
	if (!quantize_coordinate(points[i][X], tolerance, vertex.x) ||
	    !quantize_coordinate(points[i][Y], tolerance, vertex.y) ||
	    !quantize_coordinate(points[i][Z], tolerance, vertex.z))
	    return false;
	vertex.orig_ind = i;
	vertices.push_back(vertex);
    }
    std::sort(vertices.begin(), vertices.end(), hash_vertex_less);
    return true;
}

static bool
hash_vertex_coordinates_equal(const tm_hash_vert &a, const tm_hash_vert &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static void
hash_update(struct bu_data_hash_state *state, const void *data, size_t size)
{
    bu_data_hash_update(state, data, size);
}

static bool
canonical_quantized_point(tm_quantized_point &output,
	const struct bg_pca_frame *frame, const point_t point, fastf_t tolerance)
{
    vect_t delta;
    VSUB2(delta, point, frame->center);
    return quantize_coordinate(VDOT(delta, frame->xaxis), tolerance, output.x) &&
	quantize_coordinate(VDOT(delta, frame->yaxis), tolerance, output.y) &&
	quantize_coordinate(VDOT(delta, frame->zaxis), tolerance, output.z);
}

} /* namespace */

extern "C" int
bg_trimesh_diff(
	const int *f1, size_t num_f1, const point_t *p1, size_t num_p1,
        const int *f2, size_t num_f2, const point_t *p2, size_t num_p2,
	fastf_t dist_tol)
{
    if (!(dist_tol > 0.0) || !std::isfinite(dist_tol) ||
	num_f1 != num_f2 || num_p1 != num_p2 ||
	!valid_points(p1, num_p1) || !valid_points(p2, num_p2) ||
	!valid_faces(f1, num_f1, num_p1) || !valid_faces(f2, num_f2, num_p2))
	return 1;

    std::vector<size_t> ranks1;
    std::vector<size_t> ranks2;
    if (points_match_by_index(p1, p2, num_p1, dist_tol)) {
	ranks1.resize(num_p1);
	ranks2.resize(num_p2);
	std::iota(ranks1.begin(), ranks1.end(), 0u);
	std::iota(ranks2.begin(), ranks2.end(), 0u);
    } else {
	std::vector<tm_diff_vert> vertices1;
	std::vector<tm_diff_vert> vertices2;
	sorted_vertices(vertices1, p1, num_p1);
	sorted_vertices(vertices2, p2, num_p2);
	for (size_t i = 0; i < vertices1.size(); i++) {
	    if (!vertex_equal(vertices1[i], vertices2[i], dist_tol))
		return 1;
	}
	rank_by_original_index(ranks1, vertices1);
	rank_by_original_index(ranks2, vertices2);
    }
    std::vector<tm_triangle> triangles1;
    std::vector<tm_triangle> triangles2;
    if (!canonical_faces(triangles1, f1, num_f1, ranks1) ||
	!canonical_faces(triangles2, f2, num_f2, ranks2))
	return 1;
    return triangles1 == triangles2 ? 0 : 1;
}

extern "C" unsigned long long
bg_trimesh_hash(
	const int *faces, size_t faceCount, const point_t *points, size_t pointCount,
	fastf_t tolerance)
{
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) ||
	!valid_points(points, pointCount) ||
	!valid_faces(faces, faceCount, pointCount))
	return 0;

    std::vector<tm_hash_vert> vertices;
    if (!hash_vertices(vertices, points, pointCount, tolerance))
	return 0;

    struct bu_data_hash_state *state = bu_data_hash_create();
    const uint32_t version = 2;
    hash_update(state, &version, sizeof(version));
    hash_update(state, &pointCount, sizeof(pointCount));
    hash_update(state, &faceCount, sizeof(faceCount));
    for (const tm_hash_vert &vertex : vertices) {
	hash_update(state, &vertex.x, sizeof(vertex.x));
	hash_update(state, &vertex.y, sizeof(vertex.y));
	hash_update(state, &vertex.z, sizeof(vertex.z));
    }

    bool ambiguousVertices = false;
    for (size_t i = 1; i < vertices.size(); i++) {
	if (hash_vertex_coordinates_equal(vertices[i - 1], vertices[i])) {
	    ambiguousVertices = true;
	    break;
	}
    }
    const uint8_t topologyPresent = ambiguousVertices ? 0u : 1u;
    hash_update(state, &topologyPresent, sizeof(topologyPresent));
    if (!ambiguousVertices && faceCount) {
	std::vector<size_t> ranks(pointCount);
	for (size_t i = 0; i < vertices.size(); i++)
	    ranks[vertices[i].orig_ind] = i;
	std::vector<tm_triangle> triangles;
	if (!canonical_faces(triangles, faces, faceCount, ranks)) {
	    bu_data_hash_destroy(state);
	    return 0;
	}
	for (const tm_triangle &triangle : triangles)
	    hash_update(state, triangle.data(), sizeof(size_t) * triangle.size());
    }

    const unsigned long long result = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return result;
}

extern "C" int
bg_trimesh_pca_get_signature(
	struct bg_trimesh_pca_signature *signature,
	const int *faces, size_t faceCount, const point_t *points,
	size_t pointCount, fastf_t tolerance, fastf_t min_relative_axis_gap)
{
    if (!signature || !(tolerance > 0.0) || !std::isfinite(tolerance) ||
	(faceCount && (!faces || !points)) ||
	faceCount > std::numeric_limits<size_t>::max() / 3u)
	return BRLCAD_ERROR;

    struct bg_pca_frame frame;
    if (bg_pca_canonical_frame(&frame, pointCount, points,
	min_relative_axis_gap) != BRLCAD_OK)
	return BRLCAD_ERROR;

    struct bu_data_hash_state *state = bu_data_hash_create();
    const uint32_t version = 2;
    hash_update(state, &version, sizeof(version));
    hash_update(state, &pointCount, sizeof(pointCount));
    hash_update(state, &faceCount, sizeof(faceCount));
    for (size_t i = 0; i < pointCount; i++) {
	tm_quantized_point point;
	if (!canonical_quantized_point(point, &frame, points[i], tolerance)) {
	    bu_data_hash_destroy(state);
	    return BRLCAD_ERROR;
	}
	hash_update(state, &point.x, sizeof(point.x));
	hash_update(state, &point.y, sizeof(point.y));
	hash_update(state, &point.z, sizeof(point.z));
    }
    for (size_t i = 0; i < faceCount * 3u; i++) {
	if (faces[i] < 0 || static_cast<size_t>(faces[i]) >= pointCount) {
	    bu_data_hash_destroy(state);
	    return BRLCAD_ERROR;
	}
	const uint64_t vertex = static_cast<uint64_t>(faces[i]);
	hash_update(state, &vertex, sizeof(vertex));
    }

    signature->frame = frame;
    signature->hash = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return BRLCAD_OK;
}

extern "C" int
bg_trimesh_pca_equal(
	const struct bg_trimesh_pca_signature *firstSignature,
	const int *firstFaces, size_t firstFaceCount, const point_t *firstPoints,
	size_t firstPointCount,
	const struct bg_trimesh_pca_signature *secondSignature,
	const int *secondFaces, size_t secondFaceCount, const point_t *secondPoints,
	size_t secondPointCount, fastf_t tolerance)
{
    if (!firstSignature || !secondSignature || !(tolerance > 0.0) ||
	!std::isfinite(tolerance) ||
	firstFaceCount != secondFaceCount || firstPointCount != secondPointCount ||
	!valid_points(firstPoints, firstPointCount) ||
	!valid_points(secondPoints, secondPointCount) ||
	!valid_faces(firstFaces, firstFaceCount, firstPointCount) ||
	!valid_faces(secondFaces, secondFaceCount, secondPointCount))
	return 1;

    const vect_t *firstAxes[3] = {&firstSignature->frame.xaxis,
	&firstSignature->frame.yaxis, &firstSignature->frame.zaxis};
    const vect_t *secondAxes[3] = {&secondSignature->frame.xaxis,
	&secondSignature->frame.yaxis, &secondSignature->frame.zaxis};
    for (size_t i = 0; i < firstPointCount; i++) {
	vect_t firstDelta;
	vect_t secondDelta;
	VSUB2(firstDelta, firstPoints[i], firstSignature->frame.center);
	VSUB2(secondDelta, secondPoints[i], secondSignature->frame.center);
	for (size_t axis = 0; axis < 3; axis++) {
	    if (!NEAR_EQUAL(VDOT(firstDelta, *firstAxes[axis]),
		    VDOT(secondDelta, *secondAxes[axis]), tolerance))
		return 1;
	}
    }
    for (size_t i = 0; i < firstFaceCount * 3u; i++) {
	if (firstFaces[i] != secondFaces[i])
	    return 1;
    }
    return 0;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
