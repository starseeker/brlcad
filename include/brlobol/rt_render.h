/*                    R T _ R E N D E R . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/rt_render.h */

#ifndef BRLOBOL_RT_RENDER_H
#define BRLOBOL_RT_RENDER_H

#include "brlobol/defines.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>

#include <atomic>
#include <stddef.h>
#include <stdint.h>
#include <vector>

class BRLObolViewController;
struct db_i;

/** Policy for one retained librt frame.  Values are deliberately bounded by
 * the endpoint property layer; this standalone type also validates them so
 * direct library callers fail rather than silently selecting another mode. */
struct BRLOBOL_EXPORT BRLObolRtRenderSettings {
    BRLObolRtRenderSettings(void);

    unsigned int width;
    unsigned int height;
    unsigned int workers;
    unsigned int samples;
    SbColor backgroundBottom;
    SbColor backgroundTop;
};

struct BRLOBOL_EXPORT BRLObolRtRenderStatus {
    BRLObolRtRenderStatus(void);
    void clear(void);

    uint64_t geometryRevision;
    uint64_t presentationRevision;
    uint64_t raysShot;
    unsigned int preparedSources;
    unsigned int width;
    unsigned int height;
    int complete;
    int cancelled;
};

/** Optional retained-ray outputs.  Depth is normalized camera depth in
 * [0, 1], where 1 denotes no hit.  sourceIdentity is zero for background and
 * otherwise resolves through BRLObolRtRenderer::getSourceIdentity(). */
struct BRLOBOL_EXPORT BRLObolRtRenderPlanes {
    BRLObolRtRenderPlanes(void);
    void clear(void);

    std::vector<float> depth;
    std::vector<uint32_t> sourceIdentity;
};

/** One source referenced by an RT identity plane.  Identifiers are scoped to
 * the renderer's current retained snapshot and must be resolved before its
 * next synchronize() call that changes source membership. */
struct BRLOBOL_EXPORT BRLObolRtSourceIdentity {
    BRLObolRtSourceIdentity(void);
    void clear(void);

    struct db_i *database;
    SbString instanceKey;
    SbString path;
    uint32_t sourceRevision;
};

/**
 * Retained librt scene adapter for an Obol view.
 *
 * synchronize() is called on the controller/Obol owner thread.  It copies
 * source visibility and presentation state and prepares only newly changed
 * database sources.  render() subsequently touches no Coin nodes, so an
 * endpoint may execute it on a worker and discard the result when cancellation
 * is requested.  The renderer does not own an image stream or native host.
 */
class BRLOBOL_EXPORT BRLObolRtRenderer {
public:
    BRLObolRtRenderer(void);
    ~BRLObolRtRenderer(void);

    BRLObolRtRenderer(const BRLObolRtRenderer &) = delete;
    BRLObolRtRenderer &operator=(const BRLObolRtRenderer &) = delete;

    SbBool synchronize(BRLObolViewController *controller);
    SbBool render(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb,
	BRLObolRtRenderStatus *status = NULL,
	const std::atomic_bool *cancelled = NULL);

    /** Render RGB together with optional depth and source-identity planes.
     * Plane samples use a pixel-center ray, independently of RGB sampling. */
    SbBool renderWithPlanes(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes &planes,
	BRLObolRtRenderStatus *status = NULL,
	const std::atomic_bool *cancelled = NULL);

    /** Render a contiguous range of rows into a retained RGB image.  The
     * first call for a frame initializes @p rgb to the configured background;
     * subsequent ranges preserve already completed rows.  A completed range
     * is byte-identical to the same rows from render(), which lets an endpoint
     * publish bounded progressive passes without touching Coin on its worker.
     */
    SbBool renderRows(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, unsigned int firstRow,
	unsigned int rowCount, BRLObolRtRenderStatus *status = NULL,
	const std::atomic_bool *cancelled = NULL);

    /** Range form of renderWithPlanes().  The first range of a frame must
     * start at row zero so RGB and plane backgrounds are initialized. */
    SbBool renderRowsWithPlanes(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes &planes,
	unsigned int firstRow, unsigned int rowCount,
	BRLObolRtRenderStatus *status = NULL,
	const std::atomic_bool *cancelled = NULL);

    SbBool getSourceIdentity(uint32_t identifier,
	BRLObolRtSourceIdentity &identity) const;
    void clear(void);

    unsigned int getPreparedSourceCount(void) const;
    uint64_t getGeometryRevision(void) const;
    uint64_t getPresentationRevision(void) const;

private:
    SbBool renderRowsInternal(const BRLObolRtRenderSettings &settings,
	std::vector<unsigned char> &rgb, BRLObolRtRenderPlanes *planes,
	unsigned int firstRow, unsigned int rowCount,
	BRLObolRtRenderStatus *status, const std::atomic_bool *cancelled);

    struct Private;
    Private *p;
};

#endif /* BRLOBOL_RT_RENDER_H */
