/*                       P C A . C P P
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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
/** @addtogroup pca */
/** @{ */
/** @file libbg/pca.cpp
 *
 * @brief
 * Principle Component Analysis.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <algorithm>
#include <cmath>

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push /* start new diagnostic pragma */
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#elif defined(__clang__)
#  pragma clang diagnostic push /* start new diagnostic pragma */
#  pragma clang diagnostic ignored "-Wdocumentation"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#endif
#include <Eigen/Eigenvalues>
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop /* end ignoring warnings */
#elif defined(__clang__)
#  pragma clang diagnostic pop /* end ignoring warnings */
#endif


#include "vmath.h"
#include "bn/mat.h"
#include "bg/pca.h"

static int
pca_frame(struct bg_pca_frame *frame, size_t npnts, const point_t *pnts)
{
    if (!frame || npnts == 0 || !pnts)
	return BRLCAD_ERROR;

    // 1.  Find the center point
    double center[3] = {0.0, 0.0, 0.0};
    double centerCompensation[3] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < npnts; i++) {
	if (!isfinite(pnts[i][X]) || !isfinite(pnts[i][Y]) ||
	    !isfinite(pnts[i][Z]))
	    return BRLCAD_ERROR;
	for (int axis = 0; axis < 3; axis++) {
	    const double corrected = static_cast<double>(pnts[i][axis]) -
		centerCompensation[axis];
	    const double next = center[axis] + corrected;
	    centerCompensation[axis] = (next - center[axis]) - corrected;
	    center[axis] = next;
	}
    }
    for (int axis = 0; axis < 3; axis++)
	center[axis] /= static_cast<double>(npnts);

    // 2.  Accumulate the 3x3 covariance matrix.  This is mathematically
    // equivalent to the left singular vectors of the centered 3xN matrix,
    // but avoids allocating a second copy of every input point.
    // Symmetry leaves six independent values.  Scalar compensated sums avoid
    // Eigen expression work in the per-vertex hot loop.
    double covariance[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double covarianceCompensation[6] = {0.0, 0.0, 0.0,
	0.0, 0.0, 0.0};
    for (size_t i = 0; i < npnts; i++) {
	const double dx = static_cast<double>(pnts[i][X]) - center[X];
	const double dy = static_cast<double>(pnts[i][Y]) - center[Y];
	const double dz = static_cast<double>(pnts[i][Z]) - center[Z];
	const double values[6] = {dx * dx, dx * dy, dx * dz,
	    dy * dy, dy * dz, dz * dz};
	for (size_t valueIndex = 0; valueIndex < 6; valueIndex++) {
	    const double corrected = values[valueIndex] -
		covarianceCompensation[valueIndex];
	    const double next = covariance[valueIndex] + corrected;
	    covarianceCompensation[valueIndex] =
		(next - covariance[valueIndex]) - corrected;
	    covariance[valueIndex] = next;
	}
    }

    Eigen::Matrix3d covarianceMatrix;
    covarianceMatrix << covariance[0], covariance[1], covariance[2],
	covariance[1], covariance[3], covariance[4],
	covariance[2], covariance[4], covariance[5];

    // 3.  Solve the symmetric covariance matrix.  Eigenvalues are ascending;
    // expose the corresponding singular-value magnitudes in descending order.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covarianceMatrix);
    if (solver.info() != Eigen::Success)
	return BRLCAD_ERROR;

    // 4.  Extract the principal axes.  A 3x3 solver always gives three axes,
    // including for one- and two-point inputs.
    frame->center[X] = static_cast<fastf_t>(center[X]);
    frame->center[Y] = static_cast<fastf_t>(center[Y]);
    frame->center[Z] = static_cast<fastf_t>(center[Z]);
    const Eigen::Matrix3d axes = solver.eigenvectors();
    vect_t *outputAxes[3] = {&frame->xaxis, &frame->yaxis, &frame->zaxis};
    for (int i = 0; i < 3; i++) {
	const int axisIndex = 2 - i;
	(*outputAxes[i])[X] = axes(0, axisIndex);
	(*outputAxes[i])[Y] = axes(1, axisIndex);
	(*outputAxes[i])[Z] = axes(2, axisIndex);
	const double eigenvalue = std::max(0.0, solver.eigenvalues()(axisIndex));
	frame->singular_values[i] = static_cast<fastf_t>(std::sqrt(eigenvalue));
    }

    return BRLCAD_OK;
}

static int
pca_axis_sign(vect_t axis, const point_t center, size_t npnts,
	const point_t *pnts)
{
    fastf_t largest = 0.0;
    fastf_t sign = 0.0;
    const fastf_t relativeTolerance = 64.0 * DBL_EPSILON;
    for (size_t i = 0; i < npnts; i++) {
	vect_t delta;
	VSUB2(delta, pnts[i], center);
	const fastf_t projection = VDOT(delta, axis);
	const fastf_t magnitude = fabs(projection);
	if (magnitude > largest * (1.0 + relativeTolerance)) {
	    largest = magnitude;
	    sign = projection;
	    continue;
	}
	if (largest > 0.0 &&
	    fabs(magnitude - largest) <= largest * relativeTolerance &&
	    projection * sign < 0.0)
	    return BRLCAD_ERROR;
    }
    if (largest <= SMALL_FASTF || fabs(sign) <= SMALL_FASTF)
	return BRLCAD_ERROR;
    if (sign < 0.0)
	VREVERSE(axis, axis);
    return BRLCAD_OK;
}

extern "C" int
bg_pca_get_frame(struct bg_pca_frame *frame, size_t npnts, const point_t *pnts)
{
    return pca_frame(frame, npnts, pnts);
}

// Use SVD algorithm from Soderkvist to fit a plane to vertex points
extern "C" int
bg_pca(point_t *c, vect_t *xa, vect_t *ya, vect_t *za, size_t npnts,
	const point_t *pnts)
{
    if (!c || !xa || !ya || !za)
	return BRLCAD_ERROR;
    struct bg_pca_frame frame;
    if (pca_frame(&frame, npnts, pnts) != BRLCAD_OK)
	return BRLCAD_ERROR;
    VMOVE(*c, frame.center);
    VMOVE(*xa, frame.xaxis);
    VMOVE(*ya, frame.yaxis);
    VMOVE(*za, frame.zaxis);
    return BRLCAD_OK;
}

extern "C" int
bg_pca_canonical_frame(struct bg_pca_frame *frame, size_t npnts,
	const point_t *pnts, fastf_t min_relative_axis_gap)
{
    if (pca_frame(frame, npnts, pnts) != BRLCAD_OK)
	return BRLCAD_ERROR;

    const fastf_t minimumGap = min_relative_axis_gap > 0.0 ?
	min_relative_axis_gap : 1.0e-6;
    const fastf_t first = frame->singular_values[0];
    const fastf_t second = frame->singular_values[1];
    const fastf_t third = frame->singular_values[2];
    if (!isfinite(first) || !isfinite(second) || !isfinite(third) ||
	first <= SMALL_FASTF || second <= SMALL_FASTF ||
	first - second <= first * minimumGap ||
	(third > SMALL_FASTF && second - third <= second * minimumGap))
	return BRLCAD_ERROR;

    if (pca_axis_sign(frame->xaxis, frame->center, npnts, pnts) !=
	BRLCAD_OK ||
	pca_axis_sign(frame->yaxis, frame->center, npnts, pnts) !=
	BRLCAD_OK)
	return BRLCAD_ERROR;

    VCROSS(frame->zaxis, frame->xaxis, frame->yaxis);
    if (MAGNITUDE(frame->zaxis) <= SMALL_FASTF)
	return BRLCAD_ERROR;
    VUNITIZE(frame->zaxis);
    VCROSS(frame->yaxis, frame->zaxis, frame->xaxis);
    if (MAGNITUDE(frame->yaxis) <= SMALL_FASTF)
	return BRLCAD_ERROR;
    VUNITIZE(frame->yaxis);
    return BRLCAD_OK;
}

extern "C" int
bg_pca_frame_to_matrix(mat_t matrix, const struct bg_pca_frame *frame)
{
    if (!matrix || !frame || MAGNITUDE(frame->xaxis) <= SMALL_FASTF ||
	MAGNITUDE(frame->yaxis) <= SMALL_FASTF ||
	MAGNITUDE(frame->zaxis) <= SMALL_FASTF)
	return BRLCAD_ERROR;

    mat_t rotation;
    mat_t translation;
    MAT_IDN(rotation);
    VMOVE(&rotation[0], frame->xaxis);
    VMOVE(&rotation[4], frame->yaxis);
    VMOVE(&rotation[8], frame->zaxis);
    MAT_IDN(translation);
    MAT_DELTAS_VEC_NEG(translation, frame->center);
    bn_mat_mul(matrix, rotation, translation);
    return BRLCAD_OK;
}

extern "C" int
bg_pca_frame_relative_matrix(mat_t matrix,
	const struct bg_pca_frame *sourceFrame,
	const struct bg_pca_frame *targetFrame)
{
    if (!matrix || !sourceFrame || !targetFrame)
	return BRLCAD_ERROR;

    mat_t sourceToCanonical;
    mat_t targetToCanonical;
    mat_t canonicalToTarget;
    if (bg_pca_frame_to_matrix(sourceToCanonical, sourceFrame) != BRLCAD_OK ||
	bg_pca_frame_to_matrix(targetToCanonical, targetFrame) != BRLCAD_OK ||
	!bn_mat_inverse(canonicalToTarget, targetToCanonical))
	return BRLCAD_ERROR;
    bn_mat_mul(matrix, canonicalToTarget, sourceToCanonical);
    return BRLCAD_OK;
}

/** @} */
// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
