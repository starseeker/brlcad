/* L O D _ V I E W _ S N A P S H O T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_VIEW_SNAPSHOT_PRIVATE_H
#define LIBBOBOL_LOD_VIEW_SNAPSHOT_PRIVATE_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

static inline bool
bobol_lod_exact_scalar(double left, double right)
{
    if (!std::isfinite(left) || !std::isfinite(right))
	return false;
    uint64_t leftBits = 0;
    uint64_t rightBits = 0;
    std::memcpy(&leftBits, &left, sizeof(leftBits));
    std::memcpy(&rightBits, &right, sizeof(rightBits));
    if ((leftBits << 1) == 0)
	leftBits = 0;
    if ((rightBits << 1) == 0)
	rightBits = 0;
    return leftBits == rightBits;
}

static inline bool
bobol_lod_exact_scalar(float left, float right)
{
    if (!std::isfinite(left) || !std::isfinite(right))
	return false;
    uint32_t leftBits = 0;
    uint32_t rightBits = 0;
    std::memcpy(&leftBits, &left, sizeof(leftBits));
    std::memcpy(&rightBits, &right, sizeof(rightBits));
    if ((leftBits << 1) == 0)
	leftBits = 0;
    if ((rightBits << 1) == 0)
	rightBits = 0;
    return leftBits == rightBits;
}

/*
 * Exact camera and camera-scale identities used by the LoD coordinator.
 * Keep these allocation-free values beside the policies which consume them:
 * an exact-view quality proof must compare the same semantic inputs as a
 * camera epoch, rather than a lossy hash or an approximate pose.
 */
struct BObolLodViewSnapshot {
    bool same(const BObolLodViewSnapshot &other) const
    {
	/* Exact numeric equality intentionally treats +0 and -0 as the same
	 * camera while rejecting NaNs.  A byte comparison did the opposite and
	 * caused semantically identical restored views to miss history. */
	if (this->haveCamera != other.haveCamera ||
	    this->width != other.width || this->height != other.height ||
	    !bobol_lod_exact_scalar(this->size, other.size) ||
	    !bobol_lod_exact_scalar(this->lodScale, other.lodScale) ||
	    !bobol_lod_exact_scalar(this->curveScale, other.curveScale) ||
	    !bobol_lod_exact_scalar(this->pointScale, other.pointScale) ||
	    this->botThreshold != other.botThreshold)
	    return false;
	for (size_t i = 0; i < 16; ++i) {
	    if (!bobol_lod_exact_scalar(this->viewVolumeMatrix[i],
		    other.viewVolumeMatrix[i]))
		return false;
	}
	return true;
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    float viewVolumeMatrix[16] = {};
};

struct BObolLodViewScaleSnapshot {
    bool same(const BObolLodViewScaleSnapshot &other) const
    {
	return this->haveCamera == other.haveCamera &&
	    this->width == other.width &&
	    this->height == other.height &&
	    bobol_lod_exact_scalar(this->size, other.size) &&
	    bobol_lod_exact_scalar(this->lodScale, other.lodScale) &&
	    bobol_lod_exact_scalar(this->curveScale, other.curveScale) &&
	    bobol_lod_exact_scalar(this->pointScale, other.pointScale) &&
	    this->botThreshold == other.botThreshold &&
	    this->viewportWidth == other.viewportWidth &&
	    this->viewportHeight == other.viewportHeight &&
	    this->cameraTypeKey == other.cameraTypeKey &&
	    bobol_lod_exact_scalar(this->aspectRatio, other.aspectRatio) &&
	    bobol_lod_exact_scalar(this->focalDistance,
		other.focalDistance) &&
	    bobol_lod_exact_scalar(this->projectionScale,
		other.projectionScale);
    }

    uint8_t haveCamera = 0;
    int32_t width = 0;
    int32_t height = 0;
    double size = 0.0;
    double lodScale = 0.0;
    double curveScale = 0.0;
    double pointScale = 0.0;
    uint32_t botThreshold = 0;
    int16_t viewportWidth = 0;
    int16_t viewportHeight = 0;
    uint64_t cameraTypeKey = 0;
    float aspectRatio = 0.0f;
    float focalDistance = 0.0f;
    float projectionScale = 0.0f;
};

static_assert(std::is_trivially_copyable<BObolLodViewSnapshot>::value,
    "view signatures must remain allocation-free values");
static_assert(std::is_trivially_copyable<BObolLodViewScaleSnapshot>::value,
    "view scale signatures must remain allocation-free values");

#endif /* LIBBOBOL_LOD_VIEW_SNAPSHOT_PRIVATE_H */
