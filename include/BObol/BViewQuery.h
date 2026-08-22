/*                  B V I E W Q U E R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BViewQuery.h */

#ifndef BOBOL_BVIEWQUERY_H
#define BOBOL_BVIEWQUERY_H

#include "BObol/BDefines.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSnapAction.h"

#include <Inventor/SbVec3f.h>

#include <float.h>
#include <vector>

class BObolViewController;

/**
 * One raw scene hit produced by the common Obol view picker.  The record
 * deliberately retains the typed detail so callers may apply their own
 * feature, editing, or selection policy without re-traversing the scene.
 */
struct BOBOL_EXPORT BObolViewPickRecord {
    BObolViewPickRecord(void);

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
BOBOL_EXPORT int
bobol_view_pick_point(BObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records,
	int *submittedSourceRequestCount = NULL);

/**
 * Pick the currently displayed representation only.  This is the object and
 * path selection contract: an already visible occurrence must be selectable
 * immediately without requesting exact mesh data or constructing librt
 * acceleration state.  Editing and measurement should continue to use
 * bobol_view_pick_point when exact primitive identity is required.
 */
BOBOL_EXPORT int
bobol_view_pick_display_point(BObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records);

/**
 * Select displayed occurrence identities whose conservative projected bounds
 * overlap the top-left-origin viewport rectangle.  Compact CAD sources are
 * queried directly in O(occurrences), with no triangle export or ray grid.
 * Returns -1 when the scene has no compact occurrence source and the caller
 * should use its non-CAD display fallback.
 */
BOBOL_EXPORT int
bobol_view_pick_rectangle(BObolViewController *controller,
	int x0,
	int y0,
	int x1,
	int y1,
	std::vector<BObolViewPickRecord> &records);

/** Whole-source bounds fallback for non-mesh/noncompact database nodes. */
BOBOL_EXPORT int
bobol_view_pick_source_rectangle(BObolViewController *controller,
	int x0,
	int y0,
	int x1,
	int y1,
	std::vector<BObolViewPickRecord> &records);

/**
 * Pick the current Obol view using an explicit model-space ray.  See
 * bobol_view_pick_point for result and deferred-source semantics.
 */
BOBOL_EXPORT int
bobol_view_pick_ray(BObolViewController *controller,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records,
	int *submittedSourceRequestCount = NULL);

/**
 * One typed snap result produced by the common Obol view query API.  A caller
 * may keep application-specific snap state without depending on the action
 * object's traversal lifetime.
 */
struct BOBOL_EXPORT BObolViewSnapRecord {
    BObolViewSnapRecord(void);

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
BOBOL_EXPORT int
bobol_view_snap_point(BObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BObolViewSnapRecord &record);

/** As bobol_view_snap_point, restricted by SoBRLSnapAction::SourceFilter. */
BOBOL_EXPORT int
bobol_view_snap_point_filtered(BObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	uint32_t sourceFilter,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BObolViewSnapRecord &record);

#endif /* BOBOL_BVIEWQUERY_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
