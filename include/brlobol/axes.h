/*                         A X E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/axes.h */

#ifndef BRLOBOL_AXES_H
#define BRLOBOL_AXES_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;

class BRLOBOL_EXPORT SoBRLAxes : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLAxes);

public:
    SoSFString overlayId;
    SoSFVec3f origin;
    SoSFFloat size;
    SoSFBool visible;

    SoBRLAxes(void);
    static void initClass(void);

    SoBRLVListShape *rebuildGeometry(void);
    SoBRLVListShape *getGeometryShape(void) const;

protected:
    virtual ~SoBRLAxes(void);
};

#endif /* BRLOBOL_AXES_H */
