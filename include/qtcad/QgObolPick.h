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
#include <vector>

class QgView;

struct QTCAD_EXPORT QgObolPickRecord {
    enum PrimitiveKind {
	UNKNOWN = 0,
	LINE_SEGMENT = 1,
	POINT = 2,
	FACE = 3
    };

    QgObolPickRecord(void);

    std::string path;
    std::string sourceName;
    std::string sourceType;
    std::string materialShader;
    SbVec3f point;
    SbColor materialColor;
    uint32_t sourceId;
    int regionId;
    int airCode;
    int materialId;
    int los;
    int primitiveKind;
    int primitiveIndex;
    bool materialColorValid;
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
	std::vector<QgObolPickRecord> &records);

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
	std::vector<QgObolPickRecord> &records);

#endif /* QGOBOLPICK_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
