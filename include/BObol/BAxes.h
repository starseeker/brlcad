/*                       B A X E S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BAxes.h */

#ifndef BOBOL_BAXES_H
#define BOBOL_BAXES_H

#include "BObol/BDefines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;
struct bv_axes_state;

class BOBOL_EXPORT SoBRLAxes : public SoSeparator {
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

/** Publish axes view input and replacement geometry as one identity-preserving
 * transaction. Invalid input leaves the axes node unchanged. Attached callers
 * use this operation instead of assigning several fields before rebuilding. */
BOBOL_EXPORT extern int bobol_axes_configure_from_view(
    SoBRLAxes *axes,
    const struct bv_axes_state *state);

#endif /* BOBOL_BAXES_H */
