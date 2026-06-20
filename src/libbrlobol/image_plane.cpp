/*              I M A G E _ P L A N E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/image_plane.h"
#include "brlobol/image_source.h"
#include "image_display_util.h"

#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoTexture2.h>

SO_NODE_SOURCE(SoBRLImagePlane);

SoBRLImagePlane::SoBRLImagePlane(void) :
    textureNode(NULL),
    imageFaceSet(NULL)
{
    SO_NODE_CONSTRUCTOR(SoBRLImagePlane);

    SO_NODE_DEFINE_ENUM_VALUE(SizeMode, IMAGE_ASPECT);
    SO_NODE_DEFINE_ENUM_VALUE(SizeMode, EXPLICIT_SIZE);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, STRETCH);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, CONTAIN);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, COVER);

    SO_NODE_ADD_FIELD(imageSource, (NULL));
    SO_NODE_ADD_FIELD(sourcePath, (""));
    SO_NODE_ADD_FIELD(sourceId, (0));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(selectable, (TRUE));
    SO_NODE_ADD_FIELD(width, (1.0f));
    SO_NODE_ADD_FIELD(height, (1.0f));
    SO_NODE_ADD_FIELD(sizeMode, (IMAGE_ASPECT));
    SO_NODE_SET_SF_ENUM_TYPE(sizeMode, SizeMode);
    SO_NODE_ADD_FIELD(fit, (STRETCH));
    SO_NODE_SET_SF_ENUM_TYPE(fit, Fit);
    SO_NODE_ADD_FIELD(opacity, (1.0f));
    SO_NODE_ADD_FIELD(doubleSided, (TRUE));
    SO_NODE_ADD_FIELD(depthTest, (TRUE));
    SO_NODE_ADD_FIELD(depthWrite, (TRUE));
    SO_NODE_ADD_FIELD(realizedDataRevision, (0));
    SO_NODE_ADD_FIELD(realizedDirtyRevision, (0));
}

SoBRLImagePlane::~SoBRLImagePlane(void)
{
}

void
SoBRLImagePlane::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLImagePlane, SoSeparator, "Separator");
}

SoBRLImageSource *
SoBRLImagePlane::getImageSource(void) const
{
    SoNode *node = this->imageSource.getValue();
    if (!node || !node->isOfType(SoBRLImageSource::getClassTypeId()))
	return NULL;
    return static_cast<SoBRLImageSource *>(node);
}

int
SoBRLImagePlane::rebuildGeometry(void)
{
    this->textureNode = NULL;
    this->imageFaceSet = NULL;
    this->removeAllChildren();

    if (!this->visible.getValue())
	return 0;

    SoBRLImageSource *source = this->getImageSource();
    struct brlobol_image_payload payload;
    if (brlobol_image_payload_load(source, &payload) != 0)
	return -1;

    float requestedWidth = this->width.getValue();
    float requestedHeight = this->height.getValue();
    if (this->sizeMode.getValue() == IMAGE_ASPECT) {
	requestedWidth = requestedWidth > 0.0f ? requestedWidth : 1.0f;
	requestedHeight = requestedWidth * (float)payload.height / (float)payload.width;
    }

    float displayWidth = 0.0f;
    float displayHeight = 0.0f;
    brlobol_image_fit_size((float)payload.width, (float)payload.height,
	    requestedWidth, requestedHeight, this->fit.getValue() + 1,
	    this->sizeMode.getValue() == IMAGE_ASPECT,
	    &displayWidth, &displayHeight);

    float u0 = 0.0f;
    float v0 = 0.0f;
    float u1 = 1.0f;
    float v1 = 1.0f;
    brlobol_image_texture_rect((float)payload.width, (float)payload.height,
	    SbVec2f(-1.0f, -1.0f), 1.0f, this->fit.getValue() + 1,
	    &u0, &v0, &u1, &v1);

    SoTexture2 *texture = NULL;
    SoFaceSet *face = NULL;
    SoSeparator *quad = brlobol_image_make_textured_quad(&payload,
	    displayWidth * -0.5f, displayHeight * -0.5f, 0.0f,
	    displayWidth, displayHeight, u0, v0, u1, v1,
	    this->opacity.getValue(),
	    this->selectable.getValue(),
	    this->depthTest.getValue(),
	    this->depthWrite.getValue(),
	    this->doubleSided.getValue(),
	    &texture, &face);
    if (!quad)
	return -1;

    this->addChild(quad);
    this->textureNode = texture;
    this->imageFaceSet = face;
    this->realizedDataRevision = payload.dataRevision;
    this->realizedDirtyRevision = payload.dirtyRevision;
    return 0;
}

int
SoBRLImagePlane::syncFromSource(void)
{
    SoBRLImageSource *source = this->getImageSource();
    if (!source)
	return -1;
    if (source->refreshFromStream() != 0)
	return -1;
    if (this->getNumChildren() == 0 ||
	    this->realizedDataRevision.getValue() != source->dataRevision.getValue() ||
	    this->realizedDirtyRevision.getValue() != source->dirtyRevision.getValue())
	return this->rebuildGeometry();
    return 0;
}

SoTexture2 *
SoBRLImagePlane::getTextureNode(void) const
{
    return this->textureNode;
}

SoFaceSet *
SoBRLImagePlane::getImageFaceSet(void) const
{
    return this->imageFaceSet;
}
