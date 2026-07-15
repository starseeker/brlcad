/*                   V I E W _ Q U E R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/view_query.h */

#ifndef BRLOBOL_VIEW_QUERY_H
#define BRLOBOL_VIEW_QUERY_H

#include "brlobol/defines.h"
#include "brlobol/pick_detail.h"
#include "brlobol/snap_action.h"

#include <Inventor/SbVec3f.h>

#include <float.h>
#include <vector>

class BRLObolViewController;

/**
 * One raw scene hit produced by the common Obol view picker.  The record
 * deliberately retains the typed detail so callers may apply their own
 * feature, editing, or selection policy without re-traversing the scene.
 */
struct BRLOBOL_EXPORT BRLObolViewPickRecord {
    BRLObolViewPickRecord(void);

    SoBRLPickDetail detail;
    SbVec3f point;
    float distance;
};

/**
 * Pick the current Obol view at top-left viewport coordinates.  The caller
 * owns viewport and camera synchronization before invoking this routine.
 *
 * The result combines display traversal, source-backed exact mesh picks, and
 * librt exact picks.  @p submittedSourceRequestCount reports asynchronous
 * full-detail requests so callers can defer legacy selection fallbacks.
 */
BRLOBOL_EXPORT int
brlobol_view_pick_point(BRLObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records,
	int *submittedSourceRequestCount = NULL);

/**
 * Pick the current Obol view using an explicit model-space ray.  See
 * brlobol_view_pick_point for result and deferred-source semantics.
 */
BRLOBOL_EXPORT int
brlobol_view_pick_ray(BRLObolViewController *controller,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records,
	int *submittedSourceRequestCount = NULL);

/**
 * One typed snap result produced by the common Obol view query API.  A caller
 * may keep application-specific snap state without depending on the action
 * object's traversal lifetime.
 */
struct BRLOBOL_EXPORT BRLObolViewSnapRecord {
    BRLObolViewSnapRecord(void);

    SbVec3f point;
    SbString path;
    SbString editIntentId;
    SbString editIntentRole;
    SoBRLSnapAction::SnapKind kind;
    int primitiveIndex;
    int vertexIndex;
    int edgeSlot;
    int edgeVertexIndexA;
    int edgeVertexIndexB;
    int submittedSourceRequestCount;
    float distance;
    SbBool sourceFullDetailPending;
};

/**
 * Query an Obol view for one snap candidate.  A full-detail query can consume
 * ready source-backed results and submit bounded nonblocking requests when
 * @p consumeSourceFullDetail is true.  Callers retain all input and editing
 * policy; this function owns only retained-scene traversal and refinement.
 */
BRLOBOL_EXPORT int
brlobol_view_snap_point(BRLObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BRLObolViewSnapRecord &record);

#endif /* BRLOBOL_VIEW_QUERY_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
