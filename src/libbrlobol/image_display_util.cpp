/*          I M A G E _ D I S P L A Y _ U T I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "image_display_util.h"

#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPickStyle.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoTextureCoordinate2.h>

#include <algorithm>
#include <limits.h>

static uint32_t
image_display_clamp_u32(uint64_t value)
{
    return value > UINT32_MAX ? UINT32_MAX : (uint32_t)value;
}

static float
positive_or(float value, float fallback)
{
    return value > 0.0f ? value : fallback;
}

static float
image_display_clamp_float(float value, float minValue, float maxValue)
{
    if (value < minValue)
	return minValue;
    if (value > maxValue)
	return maxValue;
    return value;
}

int
brlobol_image_payload_load(SoBRLImageSource *source, struct brlobol_image_payload *payload)
{
    if (!source || !payload)
	return -1;

    if (source->refreshFromStream() != 0)
	return -1;

    imgstream_t *stream = source->getStream();
    if (!stream)
	return -1;

    struct imgstream_info info;
    if (imgstream_get_info(stream, &info) != 0)
	return -1;
    if (info.width == 0 || info.height == 0 || info.channels == 0)
	return -1;
    if (info.width > (size_t)INT_MAX || info.height > (size_t)INT_MAX ||
	info.channels > (size_t)INT_MAX)
	return -1;

    const size_t stride = info.width * info.channels;
    if (info.height > SIZE_MAX / stride)
	return -1;

    payload->pixels.assign(stride * info.height, 0);
    if (imgstream_read(stream, payload->pixels.data(), stride) != 0) {
	payload->pixels.clear();
	return -1;
    }

    payload->width = (int)info.width;
    payload->height = (int)info.height;
    payload->channels = (int)info.channels;
    payload->dataRevision = image_display_clamp_u32(info.generation);
    payload->dirtyRevision = image_display_clamp_u32(source->dirtyRevision.getValue());
    return 0;
}

void
brlobol_image_fit_size(float sourceWidth, float sourceHeight,
		       float requestedWidth, float requestedHeight,
		       int fit, bool preserveAspect,
		       float *displayWidth, float *displayHeight)
{
    const float nativeWidth = positive_or(sourceWidth, 1.0f);
    const float nativeHeight = positive_or(sourceHeight, 1.0f);
    float outWidth = positive_or(requestedWidth, nativeWidth);
    float outHeight = positive_or(requestedHeight, nativeHeight);

    if (fit == 0) {
	outWidth = nativeWidth;
	outHeight = nativeHeight;
    } else if (fit == 2 || preserveAspect) {
	float scale = std::min(outWidth / nativeWidth, outHeight / nativeHeight);
	if (scale <= 0.0f)
	    scale = 1.0f;
	outWidth = nativeWidth * scale;
	outHeight = nativeHeight * scale;
    } else if (fit == 3) {
	float scale = std::max(outWidth / nativeWidth, outHeight / nativeHeight);
	if (scale <= 0.0f)
	    scale = 1.0f;
	outWidth = nativeWidth * scale;
	outHeight = nativeHeight * scale;
    }

    if (displayWidth)
	*displayWidth = outWidth;
    if (displayHeight)
	*displayHeight = outHeight;
}

void
brlobol_image_texture_rect(float sourceWidth, float sourceHeight,
			   const SbVec2f &sourceCenter, float sourceZoom, int fit,
			   float *u0, float *v0, float *u1, float *v1)
{
    const float width = positive_or(sourceWidth, 1.0f);
    const float height = positive_or(sourceHeight, 1.0f);
    float zoom = positive_or(sourceZoom, 1.0f);
    if (zoom < 1.0f)
	zoom = 1.0f;

    float spanU = 1.0f / zoom;
    float spanV = 1.0f / zoom;
    float centerX = sourceCenter[0] >= 0.0f ? sourceCenter[0] : width * 0.5f;
    float centerY = sourceCenter[1] >= 0.0f ? sourceCenter[1] : height * 0.5f;
    float centerU = image_display_clamp_float(centerX / width, 0.0f, 1.0f);
    float centerV = image_display_clamp_float(centerY / height, 0.0f, 1.0f);

    float outU0 = centerU - spanU * 0.5f;
    float outV0 = centerV - spanV * 0.5f;
    outU0 = image_display_clamp_float(outU0, 0.0f, 1.0f - spanU);
    outV0 = image_display_clamp_float(outV0, 0.0f, 1.0f - spanV);

    if (fit == 1) {
	outU0 = 0.0f;
	outV0 = 0.0f;
	spanU = 1.0f;
	spanV = 1.0f;
    }

    if (u0)
	*u0 = outU0;
    if (v0)
	*v0 = outV0;
    if (u1)
	*u1 = outU0 + spanU;
    if (v1)
	*v1 = outV0 + spanV;
}

SoSeparator *
brlobol_image_make_textured_quad(const struct brlobol_image_payload *payload,
				 float x0, float y0, float z0, float width, float height,
				 float u0, float v0, float u1, float v1,
				 float opacity, SbBool selectable, SbBool depthTest,
				 SbBool depthWrite, SbBool doubleSided,
				 SoTexture2 **textureOut, SoFaceSet **faceOut)
{
    if (textureOut)
	*textureOut = NULL;
    if (faceOut)
	*faceOut = NULL;
    if (!payload || payload->pixels.empty() || payload->width <= 0 ||
	payload->height <= 0 || payload->channels <= 0)
	return NULL;

    SoSeparator *root = new SoSeparator;

    SoPickStyle *pick = new SoPickStyle;
    pick->style = selectable ? SoPickStyle::SHAPE : SoPickStyle::UNPICKABLE;
    root->addChild(pick);

    SoDepthBuffer *depth = new SoDepthBuffer;
    depth->test = depthTest;
    depth->write = depthWrite;
    depth->function = SoDepthBuffer::LEQUAL;
    root->addChild(depth);

    SoShapeHints *hints = new SoShapeHints;
    hints->vertexOrdering = SoShapeHints::COUNTERCLOCKWISE;
    hints->shapeType = doubleSided ? SoShapeHints::UNKNOWN_SHAPE_TYPE : SoShapeHints::SOLID;
    hints->faceType = SoShapeHints::CONVEX;
    root->addChild(hints);

    SoTexture2 *texture = new SoTexture2;
    texture->model = SoTexture2::REPLACE;
    texture->wrapS = SoTexture2::CLAMP;
    texture->wrapT = SoTexture2::CLAMP;
    texture->setImageData(payload->width, payload->height, payload->channels,
			  payload->pixels.data());
    root->addChild(texture);

    SoMaterial *material = new SoMaterial;
    material->diffuseColor.setValue(1.0f, 1.0f, 1.0f);
    material->transparency.set1Value(0, 1.0f - image_display_clamp_float(opacity, 0.0f, 1.0f));
    root->addChild(material);

    SoTextureCoordinate2 *texCoord = new SoTextureCoordinate2;
    texCoord->point.set1Value(0, SbVec2f(u0, v0));
    texCoord->point.set1Value(1, SbVec2f(u1, v0));
    texCoord->point.set1Value(2, SbVec2f(u1, v1));
    texCoord->point.set1Value(3, SbVec2f(u0, v1));
    root->addChild(texCoord);

    SoCoordinate3 *coords = new SoCoordinate3;
    coords->point.set1Value(0, SbVec3f(x0, y0, z0));
    coords->point.set1Value(1, SbVec3f(x0 + width, y0, z0));
    coords->point.set1Value(2, SbVec3f(x0 + width, y0 + height, z0));
    coords->point.set1Value(3, SbVec3f(x0, y0 + height, z0));
    root->addChild(coords);

    SoFaceSet *face = new SoFaceSet;
    face->numVertices.setValue(4);
    root->addChild(face);

    if (textureOut)
	*textureOut = texture;
    if (faceOut)
	*faceOut = face;
    return root;
}
