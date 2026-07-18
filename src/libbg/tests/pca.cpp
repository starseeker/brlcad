/*                            P C A . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file tests/pca.cpp */

#include "common.h"

#include <math.h>

#include "bu/app.h"
#include "bu/log.h"
#include "bg/pca.h"
#include "bg/trimesh.h"

static int
expect(const char *name, int condition)
{
    if (condition)
	return 0;
    bu_log("FAIL %s\n", name);
    return 1;
}

static void
frame_coordinates(fastf_t coordinates[3], const struct bg_pca_frame *frame,
	const point_t point)
{
    vect_t delta;
    VSUB2(delta, point, frame->center);
    coordinates[X] = VDOT(delta, frame->xaxis);
    coordinates[Y] = VDOT(delta, frame->yaxis);
    coordinates[Z] = VDOT(delta, frame->zaxis);
}

int
main(int UNUSED(argc), const char **argv)
{
    bu_setprogname(argv[0]);
    int failures = 0;

    point_t asymmetric[5] = {
	{0.0, 0.0, 0.0},
	{5.0, 0.0, 0.0},
	{0.0, 2.0, 0.0},
	{0.0, 0.0, 1.0},
	{1.0, 0.5, 0.25}
    };
    point_t transformed[5];
    for (size_t i = 0; i < 5; i++) {
	transformed[i][X] = -asymmetric[i][Y] + 11.0;
	transformed[i][Y] = asymmetric[i][X] - 7.0;
	transformed[i][Z] = asymmetric[i][Z] + 3.0;
    }

    int faces[9] = {0, 1, 2, 0, 3, 1, 1, 3, 2};
    struct bg_trimesh_pca_signature firstSignature;
    struct bg_trimesh_pca_signature secondSignature;
    failures += expect("first PCA signature", bg_trimesh_pca_get_signature(
	&firstSignature, faces, 3, asymmetric, 5, 1.0e-6, 1.0e-6) == BRLCAD_OK);
    failures += expect("rotated PCA signature", bg_trimesh_pca_get_signature(
	&secondSignature, faces, 3, transformed, 5, 1.0e-6, 1.0e-6) == BRLCAD_OK);
    failures += expect("transform-invariant signature",
	firstSignature.hash == secondSignature.hash);
    failures += expect("transform-invariant PCA equality", bg_trimesh_pca_equal(
	&firstSignature, faces, 3, asymmetric, 5, &secondSignature, faces, 3,
	transformed, 5, 1.0e-6) == 0);

    point_t distant[5];
    for (size_t i = 0; i < 5; i++) {
	distant[i][X] = -asymmetric[i][Y] + 1000000000.0;
	distant[i][Y] = asymmetric[i][X] - 1000000000.0;
	distant[i][Z] = asymmetric[i][Z] + 500000000.0;
    }
    struct bg_trimesh_pca_signature distantSignature;
    failures += expect("distant PCA signature", bg_trimesh_pca_get_signature(
	&distantSignature, faces, 3, distant, 5, 1.0e-6, 1.0e-6) == BRLCAD_OK);
    failures += expect("distant PCA equality", bg_trimesh_pca_equal(
	&firstSignature, faces, 3, asymmetric, 5, &distantSignature, faces, 3,
	distant, 5, 1.0e-6) == 0);

    mat_t firstToSecond;
    failures += expect("PCA relative transform", bg_pca_frame_relative_matrix(
	firstToSecond, &firstSignature.frame, &secondSignature.frame) == BRLCAD_OK);
    for (size_t i = 0; i < 5; i++) {
	point_t mapped;
	MAT4X3PNT(mapped, firstToSecond, asymmetric[i]);
	failures += expect("PCA relative transform maps vertices",
	    NEAR_EQUAL(mapped[X], transformed[i][X], 1.0e-8) &&
	    NEAR_EQUAL(mapped[Y], transformed[i][Y], 1.0e-8) &&
	    NEAR_EQUAL(mapped[Z], transformed[i][Z], 1.0e-8));
    }

    point_t changed[5];
    for (size_t i = 0; i < 5; i++)
	VMOVE(changed[i], transformed[i]);
    changed[4][Z] += 0.25;
    struct bg_trimesh_pca_signature changedSignature;
    failures += expect("changed PCA signature", bg_trimesh_pca_get_signature(
	&changedSignature, faces, 3, changed, 5, 1.0e-6, 1.0e-6) == BRLCAD_OK);
    failures += expect("changed signature differs",
	firstSignature.hash != changedSignature.hash);
    failures += expect("changed PCA equality rejected", bg_trimesh_pca_equal(
	&firstSignature, faces, 3, asymmetric, 5, &changedSignature, faces, 3,
	changed, 5, 1.0e-6) != 0);

    struct bg_pca_frame first;
    struct bg_pca_frame second;
    failures += expect("asymmetric canonical frame",
	bg_pca_canonical_frame(&first, 5, asymmetric, 1.0e-6) == BRLCAD_OK);
    failures += expect("rotated canonical frame",
	bg_pca_canonical_frame(&second, 5, transformed, 1.0e-6) == BRLCAD_OK);
    failures += expect("descending singular values",
	first.singular_values[0] >= first.singular_values[1] &&
	first.singular_values[1] >= first.singular_values[2]);
    for (size_t i = 0; i < 5; i++) {
	fastf_t firstCoordinates[3] = VINIT_ZERO;
	fastf_t secondCoordinates[3] = VINIT_ZERO;
	frame_coordinates(firstCoordinates, &first, asymmetric[i]);
	frame_coordinates(secondCoordinates, &second, transformed[i]);
	failures += expect("rotation-invariant coordinates",
	    NEAR_EQUAL(firstCoordinates[X], secondCoordinates[X], 1.0e-8) &&
	    NEAR_EQUAL(firstCoordinates[Y], secondCoordinates[Y], 1.0e-8) &&
	    NEAR_EQUAL(firstCoordinates[Z], secondCoordinates[Z], 1.0e-8));
    }

    point_t symmetric[8] = {
	{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0},
	{-1.0, 1.0, -1.0}, {1.0, 1.0, -1.0},
	{-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0},
	{-1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}
    };
    failures += expect("symmetric frame rejected",
	bg_pca_canonical_frame(&first, 8, symmetric, 1.0e-6) == BRLCAD_ERROR);

    point_t collinear[3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
	{2.0, 0.0, 0.0}};
    failures += expect("line frame rejected",
	bg_pca_canonical_frame(&first, 3, collinear, 1.0e-6) == BRLCAD_ERROR);
    failures += expect("one point pca", bg_pca_get_frame(&first, 1,
	collinear) == BRLCAD_OK);

    return failures ? 1 : 0;
}
