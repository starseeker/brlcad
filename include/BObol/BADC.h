/*                        B A D C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BADC.h */

#ifndef BOBOL_BADC_H
#define BOBOL_BADC_H

#include "BObol/BDefines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;
struct bv_adc_state;

class BOBOL_EXPORT SoBRLADC : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLADC);

public:
    SoSFString overlayId;
    SoSFVec3f center;
    SoSFFloat angleDegrees;
    SoSFFloat distance;
    SoSFFloat crosshairSize;
    SoSFFloat tickSize;
    SoSFColor lineColor;
    SoSFColor tickColor;
    SoSFInt32 lineWidth;
    SoSFBool visible;

    SoBRLADC(void);
    static void initClass(void);

    SoBRLVListShape *rebuildGeometry(void);
    SoBRLVListShape *getGeometryShape(void) const;
    SoBRLVListShape *getTickGeometryShape(void) const;

protected:
    virtual ~SoBRLADC(void);
};

/** Publish ADC view input and replacement geometry as one identity-preserving
 * transaction. Invalid input leaves the ADC node unchanged. Attached callers
 * use this operation instead of assigning several fields before rebuilding. */
BOBOL_EXPORT extern int bobol_adc_configure_from_view(
    SoBRLADC *adc,
    const struct bv_adc_state *state);

#endif /* BOBOL_BADC_H */
