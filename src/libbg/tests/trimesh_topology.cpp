/*             T R I M E S H _ T O P O L O G Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file tests/trimesh_topology.cpp
 *
 * Regression tests for triangle-mesh half-edge generation.  In particular,
 * the large case deliberately crosses the parallel sorting threshold used by
 * production meshes such as Lucy.
 */

#include "common.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <vector>

#include "bg.h"
#include "bu.h"

static bool
edge_less(const struct bg_trimesh_halfedge &left,
	  const struct bg_trimesh_halfedge &right)
{
    if (left.va != right.va)
	return left.va < right.va;
    if (left.vb != right.vb)
	return left.vb < right.vb;
    return left.flipped < right.flipped;
}

static struct bg_trimesh_halfedge
make_edge(int first, int second)
{
    struct bg_trimesh_halfedge edge;
    if (first < second) {
	edge.va = first;
	edge.vb = second;
	edge.flipped = 0;
    } else {
	edge.va = second;
	edge.vb = first;
	edge.flipped = 1;
    }
    return edge;
}

static int
verify_edge_generation(const std::vector<int> &faces)
{
    const int faceCount = static_cast<int>(faces.size() / 3);
    const size_t edgeCount = faces.size();
    struct bg_trimesh_halfedge *actual =
	bg_trimesh_generate_edge_list(faceCount,
	    const_cast<int *>(faces.data()));
    if (!actual) {
	bu_log("FAIL: bg_trimesh_generate_edge_list returned NULL\n");
	return 1;
    }

    for (size_t i = 1; i < edgeCount; ++i) {
	if (actual[i].va < actual[i - 1].va ||
	    (actual[i].va == actual[i - 1].va &&
	     actual[i].vb < actual[i - 1].vb)) {
	    bu_log("FAIL: generated half-edges are not globally ordered at %zu\n",
		i);
	    bu_free(actual, "generated half-edges");
	    return 1;
	}
    }

    std::vector<struct bg_trimesh_halfedge> expected;
    expected.reserve(edgeCount);
    for (int face = 0; face < faceCount; ++face) {
	const int *triangle = &faces[static_cast<size_t>(face) * 3];
	expected.push_back(make_edge(triangle[0], triangle[1]));
	expected.push_back(make_edge(triangle[1], triangle[2]));
	expected.push_back(make_edge(triangle[2], triangle[0]));
    }

    /* The public contract orders edge keys, not equal-key flipped records.
     * Canonicalize both equal-key groups before comparing every record. */
    std::sort(expected.begin(), expected.end(), edge_less);
    std::sort(actual, actual + edgeCount, edge_less);
    for (size_t i = 0; i < edgeCount; ++i) {
	if (actual[i].va != expected[i].va ||
	    actual[i].vb != expected[i].vb ||
	    actual[i].flipped != expected[i].flipped) {
	    bu_log("FAIL: generated half-edge differs at %zu\n", i);
	    bu_free(actual, "generated half-edges");
	    return 1;
	}
    }

    bu_free(actual, "generated half-edges");
    return 0;
}

static int
test_large_parallel_sort(void)
{
    /* 350,000 faces produce 1,050,000 half-edges, just beyond the production
     * parallel threshold.  A deterministic LCG distributes vertex IDs across
     * all active high-byte buckets and exercises the in-place partition. */
    const int faceCount = 350000;
    const uint32_t vertexLimit = 8u * 1024u * 1024u;
    std::vector<int> faces(static_cast<size_t>(faceCount) * 3);
    uint32_t state = 0x7f4a7c15u;
    for (int face = 0; face < faceCount; ++face) {
	int triangle[3];
	for (int corner = 0; corner < 3; ++corner) {
	    state = state * 1664525u + 1013904223u;
	    triangle[corner] = static_cast<int>(state % vertexLimit);
	}
	if (triangle[1] == triangle[0])
	    triangle[1] = (triangle[1] + 1) %
		static_cast<int>(vertexLimit);
	while (triangle[2] == triangle[0] ||
	       triangle[2] == triangle[1])
	    triangle[2] = (triangle[2] + 1) %
		static_cast<int>(vertexLimit);
	for (int corner = 0; corner < 3; ++corner)
	    faces[static_cast<size_t>(face) * 3 + corner] = triangle[corner];
    }

    if (verify_edge_generation(faces))
	return 1;
    bu_log("PASS: large parallel half-edge sort\n");
    return 0;
}

static int
test_signed_fallback_sort(void)
{
    const std::vector<int> faces = {
	3, -1, 2,
	7, 3, -4,
	2, 7, 3,
	-4, -1, 7
    };
    if (verify_edge_generation(faces))
	return 1;
    bu_log("PASS: signed-index fallback half-edge sort\n");
    return 0;
}

static int
test_face_count_overflow(void)
{
    int face[3] = {0, 1, 2};
    if (bg_trimesh_generate_edge_list(INT_MAX / 3 + 1, face) != NULL) {
	bu_log("FAIL: overflowing face count was accepted\n");
	return 1;
    }
    bu_log("PASS: overflowing face count rejected\n");
    return 0;
}

int
main(int UNUSED(argc), const char *argv[])
{
    bu_setprogname(argv[0]);
    int failures = 0;
    failures += test_signed_fallback_sort();
    failures += test_face_count_overflow();
    failures += test_large_parallel_sort();
    return failures ? 1 : 0;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
