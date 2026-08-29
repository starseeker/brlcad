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
#include <stdint.h>
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
pca_frame_from_moments(struct bg_pca_frame *frame, const double center[3],
	const double covariance[6])
{
    if (!frame || !center || !covariance)
	return BRLCAD_ERROR;
    for (size_t i = 0; i < 3; ++i)
	if (!isfinite(center[i]))
	    return BRLCAD_ERROR;
    for (size_t i = 0; i < 6; ++i)
	if (!isfinite(covariance[i]))
	    return BRLCAD_ERROR;

    Eigen::Matrix3d covarianceMatrix;
    covarianceMatrix << covariance[0], covariance[1], covariance[2],
	covariance[1], covariance[3], covariance[4],
	covariance[2], covariance[4], covariance[5];

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covarianceMatrix);
    if (solver.info() != Eigen::Success)
	return BRLCAD_ERROR;

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

static bool
pca_accumulator_finite(const struct bg_pca_accumulator *accumulator)
{
    if (!accumulator)
	return false;
    for (size_t axis = 0; axis < 3; ++axis)
	if (!isfinite(accumulator->mean[axis]))
	    return false;
    for (size_t value = 0; value < 6; ++value)
	if (!isfinite(accumulator->scatter[value]))
	    return false;
    return true;
}

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

    return pca_frame_from_moments(frame, center, covariance);
}

extern "C" void
bg_pca_accumulator_init(struct bg_pca_accumulator *accumulator)
{
    if (!accumulator)
	return;
    *accumulator = BG_PCA_ACCUMULATOR_INIT;
}

extern "C" int
bg_pca_accumulator_add(struct bg_pca_accumulator *accumulator,
	const point_t point)
{
    if (!accumulator || !point || accumulator->point_count == UINT64_MAX ||
	!isfinite(point[X]) || !isfinite(point[Y]) || !isfinite(point[Z]))
	return BRLCAD_ERROR;

    struct bg_pca_accumulator candidate = *accumulator;
    const uint64_t nextCount = candidate.point_count + 1;
    double delta[3] = {
	static_cast<double>(point[X]) - candidate.mean[X],
	static_cast<double>(point[Y]) - candidate.mean[Y],
	static_cast<double>(point[Z]) - candidate.mean[Z]
    };
    const double inverseCount = 1.0 / static_cast<double>(nextCount);
    for (size_t axis = 0; axis < 3; ++axis)
	candidate.mean[axis] += delta[axis] * inverseCount;
    const double updatedDelta[3] = {
	static_cast<double>(point[X]) - candidate.mean[X],
	static_cast<double>(point[Y]) - candidate.mean[Y],
	static_cast<double>(point[Z]) - candidate.mean[Z]
    };
    const double products[6] = {
	delta[X] * updatedDelta[X], delta[X] * updatedDelta[Y],
	delta[X] * updatedDelta[Z], delta[Y] * updatedDelta[Y],
	delta[Y] * updatedDelta[Z], delta[Z] * updatedDelta[Z]
    };
    for (size_t value = 0; value < 6; ++value) {
	candidate.scatter[value] += products[value];
	if (!isfinite(candidate.scatter[value]))
	    return BRLCAD_ERROR;
    }
    candidate.point_count = nextCount;
    *accumulator = candidate;
    return BRLCAD_OK;
}

extern "C" int
bg_pca_accumulator_merge(struct bg_pca_accumulator *target,
	const struct bg_pca_accumulator *source)
{
    if (!pca_accumulator_finite(target) ||
	!pca_accumulator_finite(source))
	return BRLCAD_ERROR;
    if (!source->point_count)
	return BRLCAD_OK;
    if (!target->point_count) {
	*target = *source;
	return BRLCAD_OK;
    }
    if (source->point_count > UINT64_MAX - target->point_count)
	return BRLCAD_ERROR;

    struct bg_pca_accumulator candidate = *target;
    const uint64_t combinedCount =
	candidate.point_count + source->point_count;
    const double targetWeight = static_cast<double>(candidate.point_count);
    const double sourceWeight = static_cast<double>(source->point_count);
    const double combinedWeight = static_cast<double>(combinedCount);
    const double delta[3] = {
	source->mean[X] - candidate.mean[X],
	source->mean[Y] - candidate.mean[Y],
	source->mean[Z] - candidate.mean[Z]
    };
    const double crossWeight = targetWeight * sourceWeight / combinedWeight;
    const double products[6] = {
	delta[X] * delta[X], delta[X] * delta[Y], delta[X] * delta[Z],
	delta[Y] * delta[Y], delta[Y] * delta[Z], delta[Z] * delta[Z]
    };
    for (size_t value = 0; value < 6; ++value) {
	candidate.scatter[value] += source->scatter[value] +
	    products[value] * crossWeight;
	if (!isfinite(candidate.scatter[value]))
	    return BRLCAD_ERROR;
    }
    for (size_t axis = 0; axis < 3; ++axis) {
	candidate.mean[axis] += delta[axis] * sourceWeight / combinedWeight;
	if (!isfinite(candidate.mean[axis]))
	    return BRLCAD_ERROR;
    }
    candidate.point_count = combinedCount;
    *target = candidate;
    return BRLCAD_OK;
}

extern "C" int
bg_pca_accumulator_frame(struct bg_pca_frame *frame,
	const struct bg_pca_accumulator *accumulator)
{
    if (!frame || !accumulator || !accumulator->point_count ||
	!pca_accumulator_finite(accumulator))
	return BRLCAD_ERROR;
    return pca_frame_from_moments(frame, accumulator->mean,
	accumulator->scatter);
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
