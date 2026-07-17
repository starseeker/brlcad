/*            I M A G E _ D I S P L A Y _ U T I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_IMAGE_DISPLAY_UTIL_H
#define LIBBOBOL_IMAGE_DISPLAY_UTIL_H

#include "common.h"

#include "BObol/BImageSource.h"

#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoSeparator.h>

#include <stdint.h>
#include <vector>

class SoFaceSet;
class SoTexture2;

struct bobol_image_payload {
    int width;
    int height;
    int channels;
    uint32_t dataRevision;
    uint32_t dirtyRevision;
    std::vector<unsigned char> pixels;
};

int bobol_image_payload_load(SoBRLImageSource *source, struct bobol_image_payload *payload);
void bobol_image_fit_size(float sourceWidth, float sourceHeight,
	float requestedWidth, float requestedHeight,
	int fit, bool preserveAspect,
	float *displayWidth, float *displayHeight);
void bobol_image_texture_rect(float sourceWidth, float sourceHeight,
	const SbVec2f &sourceCenter, float sourceZoom, int fit,
	float *u0, float *v0, float *u1, float *v1);
SoSeparator *bobol_image_make_textured_quad(const struct bobol_image_payload *payload,
	float x0, float y0, float z0, float width, float height,
	float u0, float v0, float u1, float v1,
	float opacity, SbBool selectable, SbBool depthTest,
	SbBool depthWrite, SbBool doubleSided,
	SoTexture2 **textureOut, SoFaceSet **faceOut);

#endif /* LIBBOBOL_IMAGE_DISPLAY_UTIL_H */
