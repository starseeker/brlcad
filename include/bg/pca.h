/*                          P C A . H
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

/*----------------------------------------------------------------------*/
/** @addtogroup bg_pca
 *
 * Principle Component Analysis
 *
 * Calculates an XYZ coordinate system such that it aligns with the largest
 * variations in the supplied data.  Intuitively, it "aligns" with the
 * shape of the points.
 *
 * To apply the coordinate system to the input points to relocate them to
 * the aligned position in model space, you can construct matrices and
 * apply them:
 *
 * @code
 * point_t center;
 * vect_t xaxis, yaxis, zaxis;
 * bg_pca(&center, &xaxis, &yaxis, &zaxis, pntcnt, pnts);
 * mat_t R, T, RT;
 * // Rotation
 * MAT_IDN(R);
 * VMOVE(&R[0], xaxis);
 * VMOVE(&R[4], yaxis);
 * VMOVE(&R[8], zaxis);
 * // Translation
 * MAT_IDN(T);
 * MAT_DELTAS_VEC_NEG(T, center);
 * // Combine
 * bn_mat_mul(RT, R, T);
 * // Apply
 * point_t p;
 * for (size_t i = 0; i < pntcnt; i++) {
 *     MAT4X3PNT(p, RT, pnts[i]);
 *     VMOVE(pnts[i], v);
 * }
 * @endcode
 */
/** @{ */
/* @file pca.h */

#ifndef BG_PCA_H
#define BG_PCA_H

#include "common.h"
#include <stdint.h>
#include "vmath.h"
#include "bg/defines.h"

__BEGIN_DECLS

/**
 * Principal-component frame for a point set.
 *
 * singular_values are descending and describe the spread along xaxis, yaxis,
 * and zaxis.  They let callers reject geometries whose PCA frame is not
 * unique enough to use as a transform-invariant identifier.
 */
struct bg_pca_frame {
    point_t center;
    vect_t xaxis;
    vect_t yaxis;
    vect_t zaxis;
    fastf_t singular_values[3];
};

/**
 * Allocation-free, mergeable PCA source moments.
 *
 * This is the numerically stable parallel form of Welford's accumulator.
 * scatter stores the symmetric centered sum in xx, xy, xz, yy, yz, zz
 * order.  Producers may update it while streaming points and merge bounded
 * worker results without retaining the source array.
 */
struct bg_pca_accumulator {
    uint64_t point_count;
    double mean[3];
    double scatter[6];
};

#define BG_PCA_ACCUMULATOR_INIT {0, {0.0, 0.0, 0.0}, \
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}

/** Initialize an empty PCA accumulator. */
BG_EXPORT extern void bg_pca_accumulator_init(
	struct bg_pca_accumulator *accumulator);

/** Add one finite point to a PCA accumulator. */
BG_EXPORT extern int bg_pca_accumulator_add(
	struct bg_pca_accumulator *accumulator, const point_t point);

/** Merge source moments into target without retaining either point stream. */
BG_EXPORT extern int bg_pca_accumulator_merge(
	struct bg_pca_accumulator *target,
	const struct bg_pca_accumulator *source);

/** Solve a principal-component frame from accumulated source moments. */
BG_EXPORT extern int bg_pca_accumulator_frame(
	struct bg_pca_frame *frame,
	const struct bg_pca_accumulator *accumulator);


/**
 * @brief
 * Perform a Principle Component Analysis on a set of points.
 *
 * Outputs are a center point and XYZ vectors for the coordinate system.
 *
 * @param[out]	c	Origin of aligned coordinate system
 * @param[out]	xaxis	Vector of X axis of aligned coordinate system (unit length)
 * @param[out]	yaxis	Vector of Y axis of aligned coordinate system (unit length)
 * @param[out]	zaxis	Vector of Z axis of aligned coordinate system (unit length)
 * @param[in]	npnts	Number of points in input pnts array
 * @param[in]	pnts	Array of points to analyze
 *
 * @return BRLCAD_OK success
 * @return BRLCAD_ERROR if inputs are invalid or calculation fails
 *
 */
BG_EXPORT extern int bg_pca(point_t *c, vect_t *xaxis, vect_t *yaxis, vect_t *zaxis, size_t npnts, const point_t *pnts);

/**
 * Calculate a PCA frame and its singular values.
 *
 * Unlike bg_pca, this API exposes the information needed to decide whether
 * the frame is sufficiently distinct for use as a geometry cache key.
 */
BG_EXPORT extern int bg_pca_get_frame(struct bg_pca_frame *frame,
	size_t npnts, const point_t *pnts);

/**
 * Calculate a deterministic right-handed PCA frame for asymmetric geometry.
 *
 * min_relative_axis_gap is the minimum fractional separation required between
 * neighboring nonzero singular values.  Values less than or equal to zero use
 * a conservative default.  The routine returns BRLCAD_ERROR for line-like or
 * axis-symmetric point sets, where PCA cannot provide a stable orientation.
 */
BG_EXPORT extern int bg_pca_canonical_frame(struct bg_pca_frame *frame,
	size_t npnts, const point_t *pnts, fastf_t min_relative_axis_gap);

/**
 * Construct the matrix mapping source points into a PCA frame's canonical
 * coordinates.  The output follows MAT4X3PNT conventions.
 */
BG_EXPORT extern int bg_pca_frame_to_matrix(mat_t matrix,
	const struct bg_pca_frame *frame);

/**
 * Construct the rigid transform mapping a point in sourceFrame's source
 * coordinates into targetFrame's source coordinates.  Both frames must
 * describe the same canonical geometry.  This is suitable for instancing a
 * representative mesh under a transformed occurrence.
 */
BG_EXPORT extern int bg_pca_frame_relative_matrix(mat_t matrix,
	const struct bg_pca_frame *sourceFrame,
	const struct bg_pca_frame *targetFrame);

__END_DECLS

#endif  /* BG_PCA_H */
/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
