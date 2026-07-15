/*            I M A G E _ D I S P L A Y _ U T I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBRLOBOL_IMAGE_DISPLAY_UTIL_H
#define LIBBRLOBOL_IMAGE_DISPLAY_UTIL_H

#include "common.h"

#include "brlobol/image_source.h"

#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoSeparator.h>

#include <stdint.h>
#include <vector>

class SoFaceSet;
class SoTexture2;

struct brlobol_image_payload {
    int width;
    int height;
    int channels;
    uint32_t dataRevision;
    uint32_t dirtyRevision;
    std::vector<unsigned char> pixels;
};

int brlobol_image_payload_load(SoBRLImageSource *source, struct brlobol_image_payload *payload);
void brlobol_image_fit_size(float sourceWidth, float sourceHeight,
	float requestedWidth, float requestedHeight,
	int fit, bool preserveAspect,
	float *displayWidth, float *displayHeight);
void brlobol_image_texture_rect(float sourceWidth, float sourceHeight,
	const SbVec2f &sourceCenter, float sourceZoom, int fit,
	float *u0, float *v0, float *u1, float *v1);
SoSeparator *brlobol_image_make_textured_quad(const struct brlobol_image_payload *payload,
	float x0, float y0, float z0, float width, float height,
	float u0, float v0, float u1, float v1,
	float opacity, SbBool selectable, SbBool depthTest,
	SbBool depthWrite, SbBool doubleSided,
	SoTexture2 **textureOut, SoFaceSet **faceOut);

#endif /* LIBBRLOBOL_IMAGE_DISPLAY_UTIL_H */
