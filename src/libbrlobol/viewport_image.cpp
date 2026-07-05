/*            V I E W P O R T _ I M A G E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/viewport_image.h"
#include "brlobol/image_source.h"
#include "image_display_util.h"

#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoTexture2.h>

SO_NODE_SOURCE(SoBRLViewportImage);

SoBRLViewportImage::SoBRLViewportImage(void) :
    textureNode(NULL),
    imageFaceSet(NULL)
{
    SO_NODE_CONSTRUCTOR(SoBRLViewportImage);

    SO_NODE_DEFINE_ENUM_VALUE(Layer, UNDERLAY);
    SO_NODE_DEFINE_ENUM_VALUE(Layer, OVERLAY);
    SO_NODE_DEFINE_ENUM_VALUE(Layer, HUD);
    SO_NODE_DEFINE_ENUM_VALUE(Anchor, LOWER_LEFT);
    SO_NODE_DEFINE_ENUM_VALUE(Anchor, LOWER_RIGHT);
    SO_NODE_DEFINE_ENUM_VALUE(Anchor, UPPER_LEFT);
    SO_NODE_DEFINE_ENUM_VALUE(Anchor, UPPER_RIGHT);
    SO_NODE_DEFINE_ENUM_VALUE(Anchor, CENTER);
    SO_NODE_DEFINE_ENUM_VALUE(Units, PIXELS);
    SO_NODE_DEFINE_ENUM_VALUE(Units, NORMALIZED);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, NATIVE);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, STRETCH);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, CONTAIN);
    SO_NODE_DEFINE_ENUM_VALUE(Fit, COVER);
    SO_NODE_DEFINE_ENUM_VALUE(CursorShape, CURSOR_NONE);
    SO_NODE_DEFINE_ENUM_VALUE(CursorShape, CURSOR_DEFAULT);
    SO_NODE_DEFINE_ENUM_VALUE(CursorShape, CURSOR_CROSSHAIR);
    SO_NODE_DEFINE_ENUM_VALUE(CursorShape, CURSOR_CUSTOM);

    SO_NODE_ADD_FIELD(imageSource, (NULL));
    SO_NODE_ADD_FIELD(overlayId, ("viewport::image"));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(layer, (OVERLAY));
    SO_NODE_SET_SF_ENUM_TYPE(layer, Layer);
    SO_NODE_ADD_FIELD(zOrder, (0));
    SO_NODE_ADD_FIELD(anchor, (LOWER_LEFT));
    SO_NODE_SET_SF_ENUM_TYPE(anchor, Anchor);
    SO_NODE_ADD_FIELD(units, (PIXELS));
    SO_NODE_SET_SF_ENUM_TYPE(units, Units);
    SO_NODE_ADD_FIELD(position, (0.0f, 0.0f));
    SO_NODE_ADD_FIELD(size, (0.0f, 0.0f));
    SO_NODE_ADD_FIELD(fit, (NATIVE));
    SO_NODE_SET_SF_ENUM_TYPE(fit, Fit);
    SO_NODE_ADD_FIELD(opacity, (1.0f));
    SO_NODE_ADD_FIELD(preserveAspect, (TRUE));
    SO_NODE_ADD_FIELD(sourceCenter, (-1.0f, -1.0f));
    SO_NODE_ADD_FIELD(sourceZoom, (1.0f));
    SO_NODE_ADD_FIELD(cursorVisible, (FALSE));
    SO_NODE_ADD_FIELD(cursorImagePosition, (0.0f, 0.0f));
    SO_NODE_ADD_FIELD(cursorShape, (CURSOR_NONE));
    SO_NODE_SET_SF_ENUM_TYPE(cursorShape, CursorShape);
    SO_NODE_ADD_FIELD(realizedDataRevision, (0));
    SO_NODE_ADD_FIELD(realizedDirtyRevision, (0));
}

SoBRLViewportImage::~SoBRLViewportImage(void)
{
}

void
SoBRLViewportImage::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLViewportImage, SoSeparator, "Separator");
}

SoBRLImageSource *
SoBRLViewportImage::getImageSource(void) const
{
    SoNode *node = this->imageSource.getValue();
    if (!node || !node->isOfType(SoBRLImageSource::getClassTypeId()))
	return NULL;
    return static_cast<SoBRLImageSource *>(node);
}

static void
viewport_anchor_origin(int anchor, float positionX, float positionY,
		       float width, float height, float *x0, float *y0)
{
    float outX = positionX;
    float outY = positionY;

    if (anchor == SoBRLViewportImage::LOWER_RIGHT ||
	anchor == SoBRLViewportImage::UPPER_RIGHT)
	outX -= width;
    if (anchor == SoBRLViewportImage::UPPER_LEFT ||
	anchor == SoBRLViewportImage::UPPER_RIGHT)
	outY -= height;
    if (anchor == SoBRLViewportImage::CENTER) {
	outX -= width * 0.5f;
	outY -= height * 0.5f;
    }

    if (x0)
	*x0 = outX;
    if (y0)
	*y0 = outY;
}

int
SoBRLViewportImage::rebuildGeometry(void)
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

    float displayWidth = 0.0f;
    float displayHeight = 0.0f;
    SbVec2f requested = this->size.getValue();
    brlobol_image_fit_size((float)payload.width, (float)payload.height,
			   requested[0], requested[1], this->fit.getValue(),
			   this->preserveAspect.getValue() == TRUE,
			   &displayWidth, &displayHeight);

    SbVec2f pos = this->position.getValue();
    float x0 = 0.0f;
    float y0 = 0.0f;
    viewport_anchor_origin(this->anchor.getValue(), pos[0], pos[1],
			   displayWidth, displayHeight, &x0, &y0);

    float u0 = 0.0f;
    float v0 = 0.0f;
    float u1 = 1.0f;
    float v1 = 1.0f;
    brlobol_image_texture_rect((float)payload.width, (float)payload.height,
			       this->sourceCenter.getValue(), this->sourceZoom.getValue(),
			       this->fit.getValue(), &u0, &v0, &u1, &v1);

    SoTexture2 *texture = NULL;
    SoFaceSet *face = NULL;
    SoSeparator *quad = brlobol_image_make_textured_quad(&payload,
			x0, y0, (float)this->zOrder.getValue(),
			displayWidth, displayHeight, u0, v0, u1, v1,
			this->opacity.getValue(), FALSE, FALSE, FALSE, TRUE,
			&texture, &face);
    if (!quad)
	return -1;

    SoHUDKit *hud = new SoHUDKit;
    hud->addWidget(quad);
    this->addChild(hud);

    this->textureNode = texture;
    this->imageFaceSet = face;
    this->realizedDataRevision = payload.dataRevision;
    this->realizedDirtyRevision = payload.dirtyRevision;
    return 0;
}

int
SoBRLViewportImage::syncFromSource(void)
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

SoHUDKit *
SoBRLViewportImage::getHUDKit(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoHUDKit::getClassTypeId()))
	    return static_cast<SoHUDKit *>(node);
    }
    return NULL;
}

SoTexture2 *
SoBRLViewportImage::getTextureNode(void) const
{
    return this->textureNode;
}

SoFaceSet *
SoBRLViewportImage::getImageFaceSet(void) const
{
    return this->imageFaceSet;
}
