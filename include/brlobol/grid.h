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

#include "bv.h"
#include "vmath.h"

#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFInt32.h>
#include <Inventor/fields/SoSFMatrix.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLVListShape;

class BRLOBOL_EXPORT SoBRLGrid : public SoSeparator
{
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLGrid);

public:
    SoSFString overlayId;
    SoSFVec3f center;
    SoSFFloat spacing;
    SoSFFloat spacingV;
    SoSFInt32 divisions;
    SoSFInt32 majorDivisionsH;
    SoSFInt32 majorDivisionsV;
    SoSFBool adaptive;
    SoSFBool visible;
    SoSFBool snapEnabled;
    SoSFBool overlayIntent;
    SoSFMatrix modelToView;
    SoSFFloat viewScale;
    SoSFFloat baseToLocal;
    SoSFFloat aspect;
    SoSFInt32 viewportWidth;
    SoSFInt32 viewportHeight;
    SoSFFloat targetPixelSpacing;
    SoSFFloat minimumPixelSpacing;
    SoSFColor minorColor;
    SoSFColor majorColor;
    SoSFColor axisColor;
    SoSFFloat effectiveSpacingH;
    SoSFFloat effectiveSpacingV;
    SoSFFloat pixelSpacingH;
    SoSFFloat pixelSpacingV;
    SoSFInt32 minorSegmentCount;
    SoSFInt32 majorSegmentCount;
    SoSFInt32 axisSegmentCount;

    SoBRLGrid(void);
    static void initClass(void);

    SoBRLVListShape *rebuildGeometry(void);
    SoBRLVListShape *getGeometryShape(void) const;
    SoBRLVListShape *getMinorShape(void) const;
    SoBRLVListShape *getMajorShape(void) const;
    SoBRLVListShape *getAxisShape(void) const;
    int getMinorSegmentCount(void) const;
    int getMajorSegmentCount(void) const;
    int getAxisSegmentCount(void) const;
    int getTotalSegmentCount(void) const;

protected:
    virtual ~SoBRLGrid(void);
};

BRLOBOL_EXPORT extern int brlobol_grid_configure_from_view(
    SoBRLGrid *grid,
    const struct bv_grid_state *state,
    const mat_t model2view,
    fastf_t view_scale,
    fastf_t base2local,
    int width,
    int height);

BRLOBOL_EXPORT extern int brlobol_grid_configure_from_view_context(
    SoBRLGrid *grid,
    const struct bv_grid_state *state,
    const void *view_ctx);

#endif /* BRLOBOL_GRID_H */
