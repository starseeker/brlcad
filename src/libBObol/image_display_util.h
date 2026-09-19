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
#include <Inventor/tools/SbModernUtils.h>
#include <Inventor/nodes/SoSeparator.h>

#include <stdint.h>
#include <exception>
#include <memory>
#include <vector>

class SoBRLViewportImage;
class BObolViewController;
class SoFaceSet;
class SoGroup;
class SoSFUInt32;
class SoTexture2;
struct bobol_image_payload;

enum class BObolViewportGeometryPublication {
    Preserve,
    Rebuild,
    RebuildFromRetainedTexture
};

/* Prepare a viewport's scalar presentation state and, when requested,
 * successor geometry without changing the live node. Callers may prepare
 * related parent lists, commit every participant, and then notify. */
class BObolPreparedViewportImage {
public:
    explicit BObolPreparedViewportImage(SoBRLViewportImage &target);
    BObolPreparedViewportImage(SoBRLViewportImage &target,
	int layer, SbBool visible,
	BObolViewportGeometryPublication geometryPublication =
	    BObolViewportGeometryPublication::Rebuild);
    ~BObolPreparedViewportImage();

    BObolPreparedViewportImage(const BObolPreparedViewportImage &) = delete;
    BObolPreparedViewportImage &operator=(const BObolPreparedViewportImage &) = delete;

    SoBRLViewportImage &next();
    bool prepare(BObolViewportGeometryPublication geometryPublication);
    bool prepare(const struct bobol_image_payload &payload);
    bool valid() const;
    void commit();
    void restore();
    void notify(std::exception_ptr &failure);

private:
    bool preparePayload(const struct bobol_image_payload *payload,
	bool publishGeometry);

    struct Impl;
    std::unique_ptr<Impl> impl;
};

enum class BObolFramebufferRootInsertion { First, Last };

class BObolPreparedFramebufferRoots {
public:
    BObolPreparedFramebufferRoots(BObolViewController &controller,
	SoNode *viewport, SoGroup *destination,
	BObolFramebufferRootInsertion insertion);
    ~BObolPreparedFramebufferRoots();

    BObolPreparedFramebufferRoots(const BObolPreparedFramebufferRoots &) = delete;
    BObolPreparedFramebufferRoots &operator=(const BObolPreparedFramebufferRoots &) = delete;

    void commit();
    void notify(std::exception_ptr &failure);

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

struct bobol_image_payload {
    int width;
    int height;
    int channels;
    uint32_t dataRevision;
    uint32_t dirtyRevision;
    std::vector<unsigned char> pixels;
};

int bobol_image_payload_load(SoBRLImageSource *source, struct bobol_image_payload *payload);
int bobol_image_payload_load_current(SoBRLImageSource *source,
	struct bobol_image_payload *payload);
int bobol_image_payload_load_info(imgstream_t *stream,
	const struct imgstream_info &info, uint32_t dirtyRevision,
	struct bobol_image_payload *payload);
int bobol_image_payload_load_retained(const SoBRLViewportImage *viewport,
	struct bobol_image_payload *payload);
void bobol_image_fit_size(float sourceWidth, float sourceHeight,
	float requestedWidth, float requestedHeight,
	int fit, bool preserveAspect,
	float *displayWidth, float *displayHeight);
void bobol_image_texture_rect(float sourceWidth, float sourceHeight,
	const SbVec2f &sourceCenter, float sourceZoom,
	float *u0, float *v0, float *u1, float *v1);
SbModernUtils::SoNodeRef bobol_image_make_textured_quad(const struct bobol_image_payload *payload,
	float x0, float y0, float z0, float width, float height,
	float u0, float v0, float u1, float v1,
	float opacity, SbBool selectable, SbBool depthTest,
	SbBool depthWrite, SbBool doubleSided,
	SoTexture2 **textureOut, SoFaceSet **faceOut);
SbModernUtils::SoNodeRef bobol_viewport_image_make_geometry(
	const SoBRLViewportImage &viewport,
	const struct bobol_image_payload &payload,
	SoTexture2 **textureOut, SoFaceSet **faceOut);
void bobol_image_publish_geometry(SoSeparator &target, SoNode *child,
	SoTexture2 *texture, SoFaceSet *face,
	SoTexture2 *&textureSlot, SoFaceSet *&faceSlot,
	SoSFUInt32 &dataRevision, SoSFUInt32 &dirtyRevision,
	uint32_t nextDataRevision, uint32_t nextDirtyRevision,
	bool publishRevisions);

#endif /* LIBBOBOL_IMAGE_DISPLAY_UTIL_H */
