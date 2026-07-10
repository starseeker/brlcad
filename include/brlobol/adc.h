/*                         A D C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/adc.h */

#ifndef BRLOBOL_ADC_H
#define BRLOBOL_ADC_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;

class BRLOBOL_EXPORT SoBRLADC : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLADC);

public:
    SoSFString overlayId;
    SoSFVec3f center;
    SoSFFloat angleDegrees;
    SoSFFloat distance;
    SoSFFloat crosshairSize;
    SoSFFloat tickSize;
    SoSFBool visible;

    SoBRLADC(void);
    static void initClass(void);

    SoBRLVListShape *rebuildGeometry(void);
    SoBRLVListShape *getGeometryShape(void) const;

protected:
    virtual ~SoBRLADC(void);
};

#endif /* BRLOBOL_ADC_H */
