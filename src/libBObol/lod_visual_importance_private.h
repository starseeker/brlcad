/*        L O D _ V I S U A L _ I M P O R T A N C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_VISUAL_IMPORTANCE_PRIVATE_H
#define LIBBOBOL_LOD_VISUAL_IMPORTANCE_PRIVATE_H

#include "BObol/BViewLod.h"

#include <algorithm>
#include <cmath>
#include <limits>

/* Visual priority is a policy class, not a numeric multiplier.  Keeping user
 * intent lexicographic prevents a sufficiently large ordinary object from
 * outranking the sole selection which must remain suitable for editing. */
enum BObolLodVisualEmphasis {
    BOBOL_LOD_VISUAL_ORDINARY = 0,
    BOBOL_LOD_VISUAL_HIGHLIGHTED = 1,
    BOBOL_LOD_VISUAL_SELECTED = 2
};

static constexpr double BObolLodHighlightedMaximumErrorPixels = 1.5;

inline unsigned int
bobol_lod_visual_emphasis(int emphasis)
{
    if (emphasis >= BOBOL_LOD_VISUAL_SELECTED)
	return BOBOL_LOD_VISUAL_SELECTED;
    if (emphasis == BOBOL_LOD_VISUAL_HIGHLIGHTED)
	return BOBOL_LOD_VISUAL_HIGHLIGHTED;
    return BOBOL_LOD_VISUAL_ORDINARY;
}

/* A projected box is often sparse.  Its square-root area captures occupied
 * screen scale, perimeter protects thin silhouettes, and one quarter of the
 * full diameter preserves a conservative signal for degenerate projections.
 * Projection owns viewport clipping; this function only combines its three
 * renderer-neutral measurements. */
inline double
bobol_lod_visual_footprint(double area, double perimeter, double diameter)
{
    const double safeArea = std::isfinite(area) ? std::max(0.0, area) : 0.0;
    const double safePerimeter = std::isfinite(perimeter) ?
	std::max(0.0, perimeter) : 0.0;
    const double safeDiameter = std::isfinite(diameter) ?
	std::max(0.0, diameter) : 0.0;
    const double visibleSpan = std::max(std::sqrt(safeArea),
	safePerimeter * 0.25);
    return visibleSpan > 0.0 ? visibleSpan : safeDiameter * 0.25;
}

inline double
bobol_lod_normalized_visual_error(double projectedError,
	double targetPixelError)
{
    if (!std::isfinite(projectedError) || projectedError < 0.0)
	return std::numeric_limits<double>::infinity();
    const double target = std::isfinite(targetPixelError) &&
	targetPixelError > 0.0 ? targetPixelError :
	static_cast<double>(std::numeric_limits<float>::min());
    return projectedError / target;
}

inline bool
bobol_lod_visual_prominent(double footprint)
{
    return std::isfinite(footprint) &&
	footprint >= BObolViewLodState::ProminentFootprintPixels;
}

inline double
bobol_lod_protected_visual_error(unsigned int emphasis, double footprint,
	double targetPixelError)
{
    const double target = std::isfinite(targetPixelError) &&
	targetPixelError > 0.0 ? targetPixelError :
	static_cast<double>(std::numeric_limits<float>::min());
    if (emphasis >= BOBOL_LOD_VISUAL_SELECTED)
	return target;
    if (emphasis == BOBOL_LOD_VISUAL_HIGHLIGHTED)
	return std::max(target, BObolLodHighlightedMaximumErrorPixels);
    return bobol_lod_visual_prominent(footprint) ?
	BObolViewLodState::prominentMaximumProjectedError(target) :
	std::numeric_limits<double>::infinity();
}

/* Once explicit priority and recognizability floors are satisfied, spend the
 * remaining scene budget as a rate-distortion problem.  Screen-space error
 * already includes projection scale, so footprint belongs in the affected
 * image span exactly once rather than in both error and benefit. */
inline double
bobol_lod_marginal_visual_benefit(double footprint,
	double currentProjectedError, double nextProjectedError,
	double targetPixelError)
{
    const double current = bobol_lod_normalized_visual_error(
	currentProjectedError, targetPixelError);
    const double next = bobol_lod_normalized_visual_error(
	nextProjectedError, targetPixelError);
    if (!std::isfinite(current) || !std::isfinite(next) || next >= current)
	return 0.0;
    return std::max(1.0, footprint) * (current - next);
}

#endif /* LIBBOBOL_LOD_VISUAL_IMPORTANCE_PRIVATE_H */
