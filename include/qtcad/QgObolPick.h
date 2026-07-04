/*                    Q G O B O L P I C K . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolPick.h */

#ifndef QGOBOLPICK_H
#define QGOBOLPICK_H

#include "qtcad/defines.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbVec3f.h>

#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

class QgView;

struct QTCAD_EXPORT QgObolPickRecord {
    enum PrimitiveKind {
	UNKNOWN = 0,
	LINE_SEGMENT = 1,
	POINT = 2,
	FACE = 3,
	IMPLICIT_SOLID = 4
    };

    QgObolPickRecord(void);

    std::string path;
    std::string sourceName;
    std::string sourceType;
    std::string featureName;
    std::string materialShader;
    std::string editIntentId;
    std::string editIntentRole;
    SbVec3f point;
    SbVec3f modelPoint;
    SbColor materialColor;
    float distance;
    uint32_t sourceId;
    int regionId;
    int airCode;
    int materialId;
    int los;
    int primitiveKind;
    int primitiveIndex;
    int featurePrimitiveIndex;
    int faceVertexIndexA;
    int faceVertexIndexB;
    int faceVertexIndexC;
    int nearestFaceEdgeSlot;
    int nearestFaceEdgeVertexIndexA;
    int nearestFaceEdgeVertexIndexB;
    int nearestFaceVertexSlot;
    int nearestFaceVertexIndex;
    bool materialColorValid;
    bool featurePickResolved;
    std::vector<std::pair<std::string, std::string> > featurePrimitiveMetadata;
};

/**
 * Pick BRL-CAD Obol geometry in @p display at Qt widget pixel coordinate
 * (@p x, @p y).  Qt coordinates are top-left based; this helper converts them
 * to the lower-left viewport coordinates expected by SoRayPickAction.
 *
 * Returns the number of BRL-CAD Obol pick records appended to @p records.
 */
QTCAD_EXPORT int qg_obol_pick_point(QgView *display,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount = 0);

/**
 * Pick BRL-CAD Obol geometry using an explicit model-space ray.  This helper
 * uses the same display-level, source-backed mesh, and cache-backed librt
 * exact picking paths as point picks, but lets ray-selection tools provide
 * their already computed ray directly.
 *
 * Returns the number of BRL-CAD Obol pick records appended to @p records.
 */
QTCAD_EXPORT int qg_obol_pick_ray(QgView *display,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount = 0);

/**
 * Pick BRL-CAD Obol geometry intersecting the Qt widget pixel rectangle
 * bounded by (@p x0, @p y0) and (@p x1, @p y1).  This helper samples the
 * rectangle through Obol point picks and returns unique BRL-CAD identities.
 *
 * Returns the number of BRL-CAD Obol pick records appended to @p records.
 */
QTCAD_EXPORT int qg_obol_pick_rect(QgView *display,
	int x0,
	int y0,
	int x1,
	int y1,
	float radiusPixels,
	bool firstOnly,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount = 0);

/**
 * Apply Obol feature primitive state for a resolved pick record.  Database
 * object picks normally do not resolve to feature-store primitives; this helper
 * is intended for command-result and overlay feature picks.
 *
 * Returns 1 when the record resolved to an Obol feature primitive and all
 * requested state was applied.
 */
QTCAD_EXPORT int qg_obol_pick_apply_feature_state(QgView *display,
	const QgObolPickRecord &record,
	bool select,
	bool highlight);

/**
 * Apply Obol feature primitive state for all resolved pick records, grouping
 * primitives by parent feature before replacing selected/highlighted state.
 *
 * Returns the number of Obol parent features updated.
 */
QTCAD_EXPORT int qg_obol_pick_apply_feature_states(QgView *display,
	const std::vector<QgObolPickRecord> &records,
	bool select,
	bool highlight);

#endif /* QGOBOLPICK_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
