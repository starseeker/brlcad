/*                  Q G O B O L S N A P . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolSnap.h */

#ifndef QGOBOLSNAP_H
#define QGOBOLSNAP_H

#include "qtcad/defines.h"

#include <Inventor/SbVec3f.h>

#include <stdint.h>
#include <string>

class QgView;

struct QTCAD_EXPORT QgObolSnapRecord {
    enum SnapKind {
	NONE = 0,
	ENDPOINT = 1,
	MIDPOINT = 2,
	LINE_NEAREST = 4,
	FACE_NEAREST = 8,
	CENTER = 16,
	CONSTRUCTION_PLANE = 32,
	ALL_KINDS = ENDPOINT | MIDPOINT | LINE_NEAREST | FACE_NEAREST | CENTER | CONSTRUCTION_PLANE
    };

    QgObolSnapRecord(void);

    SbVec3f point;
    std::string path;
    int kind;
    int primitiveIndex;
    float distance;
};

/**
 * Snap @p query, expressed in model/world coordinates, to Obol geometry in
 * @p display.  @p tolerance is a model-space distance.  @p enabledKinds uses
 * QgObolSnapRecord::SnapKind bits; passing 0 enables every kind.
 *
 * Returns non-zero when a snap candidate is found and @p record is populated.
 */
QTCAD_EXPORT int qg_obol_snap_point(QgView *display,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	QgObolSnapRecord &record);

#endif /* QGOBOLSNAP_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
