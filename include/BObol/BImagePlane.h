/*                 B I M A G E P L A N E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BImagePlane.h */

#ifndef BOBOL_BIMAGEPLANE_H
#define BOBOL_BIMAGEPLANE_H

#include "BObol/BDefines.h"

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/fields/SoSFUInt32.h>
#include <Inventor/nodes/SoSeparator.h>

class SoBRLImageSource;
class SoFaceSet;
class SoTexture2;

class BOBOL_EXPORT SoBRLImagePlane : public SoSeparator {
    typedef SoSeparator inherited;

    SO_NODE_HEADER(SoBRLImagePlane);

public:
    enum SizeMode {
	IMAGE_ASPECT = 0,
	EXPLICIT_SIZE = 1
    };

    enum Fit {
	STRETCH = 0,
	CONTAIN = 1,
	COVER = 2
    };

    SoSFNode imageSource;
    SoSFString sourcePath;
    SoSFUInt32 sourceId;
    SoSFBool visible;
    SoSFBool selectable;
    SoSFFloat width;
    SoSFFloat height;
    SoSFEnum sizeMode;
    SoSFEnum fit;
    SoSFFloat opacity;
    SoSFBool doubleSided;
    SoSFBool depthTest;
    SoSFBool depthWrite;
    SoSFUInt32 realizedDataRevision;
    SoSFUInt32 realizedDirtyRevision;

    SoBRLImagePlane(void);
    static void initClass(void);

    SoBRLImageSource *getImageSource(void) const;
    int rebuildGeometry(void);
    int syncFromSource(void);
    SoTexture2 *getTextureNode(void) const;
    SoFaceSet *getImageFaceSet(void) const;

protected:
    virtual ~SoBRLImagePlane(void);

private:
    SoTexture2 *textureNode;
    SoFaceSet *imageFaceSet;
};

#endif /* BOBOL_BIMAGEPLANE_H */
