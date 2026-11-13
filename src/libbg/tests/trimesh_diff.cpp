/*                 T R I M E S H _ D I F F . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file tests/trimesh_diff.cpp */

#include "common.h"

#include "bu/app.h"
#include "bu/log.h"
#include "bg/trimesh.h"

static const fastf_t test_tolerance = 1.0e-6;

static point_t base_points[4] = {
    {0.0, 0.0, 0.0},
    {2.0, 0.0, 0.0},
    {0.0, 3.0, 0.0},
    {0.0, 0.0, 4.0}
};

static int base_faces[12] = {
    0, 2, 1,
    0, 1, 3,
    1, 2, 3,
    2, 0, 3
};

static point_t permuted_points[4] = {
    {0.0, 3.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 4.0},
    {2.0, 0.0, 0.0}
};

static int permuted_faces[12] = {
    0, 2, 3,
    0, 3, 1,
    1, 2, 0,
    3, 2, 1
};

static int topology_changed_faces[12] = {
    0, 1, 3,
    0, 1, 3,
    1, 2, 3,
    2, 0, 3
};

static int reversed_winding_faces[12] = {
    0, 1, 2,
    0, 1, 3,
    1, 2, 3,
    2, 0, 3
};

static int invalid_faces[12] = {
    0, 2, 1,
    0, 1, 3,
    1, 2, 3,
    2, 0, 4
};

static int
expect(const char *name, int condition)
{
    if (condition)
	return 0;
    bu_log("FAIL %s\n", name);
    return 1;
}

int
main(int UNUSED(argc), const char **argv)
{
    bu_setprogname(argv[0]);
    int failures = 0;

    const unsigned long long baseHash = bg_trimesh_hash(base_faces, 4,
	base_points, 4, test_tolerance);
    const unsigned long long permutedHash = bg_trimesh_hash(permuted_faces,
	4, permuted_points, 4, test_tolerance);
    failures += expect("base hash", baseHash != 0);
    failures += expect("permuted equivalent mesh", bg_trimesh_diff(
	base_faces, 4, base_points, 4, permuted_faces, 4, permuted_points, 4,
	test_tolerance) == 0);
    failures += expect("permuted equivalent hash", baseHash == permutedHash);

    failures += expect("changed topology", bg_trimesh_diff(base_faces, 4,
	base_points, 4, topology_changed_faces, 4, base_points, 4,
	test_tolerance) != 0);
    failures += expect("reversed winding", bg_trimesh_diff(base_faces, 4,
	base_points, 4, reversed_winding_faces, 4, base_points, 4,
	test_tolerance) != 0);
    failures += expect("changed topology hash", baseHash != bg_trimesh_hash(
	topology_changed_faces, 4, base_points, 4, test_tolerance));

    point_t withinTolerance[4];
    point_t outsideTolerance[4];
    for (size_t i = 0; i < 4; i++) {
	VMOVE(withinTolerance[i], base_points[i]);
	VMOVE(outsideTolerance[i], base_points[i]);
    }
    withinTolerance[3][X] += test_tolerance * 0.5;
    outsideTolerance[3][X] += test_tolerance * 2.0;
    failures += expect("within tolerance", bg_trimesh_diff(base_faces, 4,
	base_points, 4, base_faces, 4, withinTolerance, 4, test_tolerance) == 0);
    failures += expect("outside tolerance", bg_trimesh_diff(base_faces, 4,
	base_points, 4, base_faces, 4, outsideTolerance, 4, test_tolerance) != 0);

    failures += expect("invalid face rejected", bg_trimesh_diff(base_faces,
	4, base_points, 4, invalid_faces, 4, base_points, 4,
	test_tolerance) != 0);
    failures += expect("invalid face hash", bg_trimesh_hash(invalid_faces,
	4, base_points, 4, test_tolerance) == 0);
    failures += expect("point-only mesh", bg_trimesh_diff(NULL, 0,
	base_points, 4, NULL, 0, permuted_points, 4, test_tolerance) == 0);
    failures += expect("point-only hash", bg_trimesh_hash(NULL, 0,
	base_points, 4, test_tolerance) != 0);

    return failures ? 1 : 0;
}
