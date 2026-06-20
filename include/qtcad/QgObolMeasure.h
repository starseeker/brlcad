/*                 Q G O B O L M E A S U R E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolMeasure.h */

#ifndef QGOBOLMEASURE_H
#define QGOBOLMEASURE_H

#include "qtcad/defines.h"

#include <Inventor/SbVec3f.h>

#include <string>

class QgView;
struct bu_color;

struct QTCAD_EXPORT QgObolMeasureGeometryRecord {
    enum PrimitiveKind {
	NONE = 0,
	LINE_SEGMENT = 1,
	FACE = 3
    };

    QgObolMeasureGeometryRecord(void);

    int shapeCount;
    int segmentCount;
    int triangleCount;
    float surfaceArea;
    float totalLength;
    SbVec3f boundsMin;
    SbVec3f boundsMax;
    bool boundsValid;
    bool hasNearestPrimitive;
    int nearestPrimitiveKind;
    int nearestPrimitiveIndex;
    std::string nearestPath;
    SbVec3f nearestPoint;
    float nearestDistance;
};

QTCAD_EXPORT int qg_obol_measure_pick_point(QgView *display,
	int x,
	int y,
	SbVec3f &point,
	std::string *path = 0);

QTCAD_EXPORT int qg_obol_measure_geometry_full_detail(QgView *display,
	const SbVec3f *query,
	QgObolMeasureGeometryRecord &record);

QTCAD_EXPORT int qg_obol_measure_update_overlay(QgView *display,
	const char *overlayId,
	const SbVec3f *points,
	int pointCount,
	const struct bu_color *color);

QTCAD_EXPORT int qg_obol_measure_clear_overlay(QgView *display,
	const char *overlayId);

#endif /* QGOBOLMEASURE_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
