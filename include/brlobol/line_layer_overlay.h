/*                 L I N E _ L A Y E R _ O V E R L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/line_layer_overlay.h */

#ifndef BRLOBOL_LINE_LAYER_OVERLAY_H
#define BRLOBOL_LINE_LAYER_OVERLAY_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

#include <stddef.h>

class SoBRLVListShape;
struct bg_line_layer_builder;

class BRLOBOL_EXPORT SoBRLLineLayerOverlay : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLLineLayerOverlay);

public:
    SoSFString overlayId;
    SoSFUInt32 sourceId;
    SoSFBool visible;
    SoSFBool selectable;

    SoBRLLineLayerOverlay(void);
    static void initClass(void);

    int rebuildGeometry(const struct bg_line_layer_builder *builder);
    SoBRLVListShape *getLayerShape(int index) const;
    int getLayerShapeCount(void) const;
    size_t getPointCount(void) const;

protected:
    virtual ~SoBRLLineLayerOverlay(void);
};

#endif /* BRLOBOL_LINE_LAYER_OVERLAY_H */
