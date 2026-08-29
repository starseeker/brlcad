/*           B L I N E L A Y E R O V E R L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLineLayerOverlay.h */

#ifndef BOBOL_BLINELAYEROVERLAY_H
#define BOBOL_BLINELAYEROVERLAY_H

#include "BObol/BDefines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

#include <stddef.h>

class SoBRLVListShape;
struct bg_line_layer_builder;

class BOBOL_EXPORT SoBRLLineLayerOverlay : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLLineLayerOverlay);

public:
    SoSFString overlayId;
    SoSFUInt32 sourceId;
    SoSFBool visible;
    SoSFBool selectable;
    SoSFBool depthTest;

    SoBRLLineLayerOverlay(void);
    static void initClass(void);

    int rebuildGeometry(const struct bg_line_layer_builder *builder);
    SoBRLVListShape *getLayerShape(int index) const;
    int getLayerShapeCount(void) const;
    size_t getPointCount(void) const;

protected:
    virtual ~SoBRLLineLayerOverlay(void);
};

#endif /* BOBOL_BLINELAYEROVERLAY_H */
