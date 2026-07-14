/*                V I E W P O R T _ I M A G E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/viewport_image.h */

#ifndef BRLOBOL_VIEWPORT_IMAGE_H
#define BRLOBOL_VIEWPORT_IMAGE_H

#include "brlobol/defines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/fields/SoSFVec2f.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLImageSource;
class SoFaceSet;
class SoHUDKit;
class SoTexture2;

class BRLOBOL_EXPORT SoBRLViewportImage : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLViewportImage);

public:
    enum Layer {
	UNDERLAY = 0,
	OVERLAY = 1,
	HUD = 2,
	INTERLAY = 3
    };

    enum Anchor {
	LOWER_LEFT = 0,
	LOWER_RIGHT = 1,
	UPPER_LEFT = 2,
	UPPER_RIGHT = 3,
	CENTER = 4
    };

    enum Units {
	PIXELS = 0,
	NORMALIZED = 1
    };

    enum Fit {
	NATIVE = 0,
	STRETCH = 1,
	CONTAIN = 2,
	COVER = 3
    };

    enum CursorShape {
	CURSOR_NONE = 0,
	CURSOR_DEFAULT = 1,
	CURSOR_CROSSHAIR = 2,
	CURSOR_CUSTOM = 3
    };

    SoSFNode imageSource;
    SoSFString overlayId;
    SoSFBool visible;
    SoSFEnum layer;
    SoSFUInt32 zOrder;
    SoSFEnum anchor;
    SoSFEnum units;
    SoSFVec2f position;
    SoSFVec2f size;
    SoSFEnum fit;
    SoSFFloat opacity;
    SoSFBool preserveAspect;
    SoSFVec2f sourceCenter;
    SoSFFloat sourceZoom;
    SoSFBool cursorVisible;
    SoSFVec2f cursorImagePosition;
    SoSFEnum cursorShape;
    SoSFUInt32 realizedDataRevision;
    SoSFUInt32 realizedDirtyRevision;

    SoBRLViewportImage(void);
    static void initClass(void);

    SoBRLImageSource *getImageSource(void) const;
    int rebuildGeometry(void);
    int syncFromSource(void);
    SoHUDKit *getHUDKit(void) const;
    SoTexture2 *getTextureNode(void) const;
    SoFaceSet *getImageFaceSet(void) const;

protected:
    virtual ~SoBRLViewportImage(void);

private:
    SoTexture2 *textureNode;
    SoFaceSet *imageFaceSet;
};

#endif /* BRLOBOL_VIEWPORT_IMAGE_H */
