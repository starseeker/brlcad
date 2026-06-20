/*                         G R I D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/grid.h */

#ifndef BRLOBOL_GRID_H
#define BRLOBOL_GRID_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;

class BRLOBOL_EXPORT SoBRLGrid : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLGrid);

public:
    SoSFString overlayId;
    SoSFVec3f center;
    SoSFFloat spacing;
    SoSFInt32 divisions;
    SoSFBool visible;

    SoBRLGrid(void);
    static void initClass(void);

    SoBRLVListShape *rebuildGeometry(void);
    SoBRLVListShape *getGeometryShape(void) const;

protected:
    virtual ~SoBRLGrid(void);
};

#endif /* BRLOBOL_GRID_H */
