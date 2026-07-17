/*              T E S T _ I M A G E _ D I S P L A Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "BObol.h"

#include "bu/log.h"
#include "imgstream/stream.h"

#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoTexture2.h>

#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

static int
texture_matches(SoTexture2 *texture, int width, int height, int channels)
{
    if (!texture)
	return 0;
    int texWidth = 0;
    int texHeight = 0;
    int texChannels = 0;
    const unsigned char *pixels = texture->getImageData(texWidth, texHeight, texChannels);
    return pixels && texWidth == width && texHeight == height && texChannels == channels;
}

static imgstream_t *
make_rgba_stream(void)
{
    imgstream_t *stream = imgstream_create(4, 2, IMGSTREAM_PIXEL_RGBA8);
    if (!stream)
	return NULL;

    unsigned char pixels[4 * 2 * 4];
    for (size_t i = 0; i < sizeof(pixels); i++)
	pixels[i] = (unsigned char)(i + 1);
    if (imgstream_write(stream, pixels, 4 * 4) != 0) {
	imgstream_destroy(stream);
	return NULL;
    }
    return stream;
}

static int
test_viewport_image(void)
{
    imgstream_t *stream = make_rgba_stream();
    CHECK(stream != NULL, "created RGBA stream for viewport image");

    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();
    CHECK(source->setStream(stream) == 0, "attached viewport source stream");

    SoBRLViewportImage *viewport = new SoBRLViewportImage;
    viewport->ref();
    viewport->imageSource.setValue(source);
    viewport->fit = SoBRLViewportImage::STRETCH;
    viewport->size.setValue(8.0f, 4.0f);
    viewport->position.setValue(10.0f, 20.0f);
    viewport->anchor = SoBRLViewportImage::LOWER_LEFT;
    viewport->sourceZoom = 1.0f;
    viewport->opacity = 0.75f;

    CHECK(viewport->rebuildGeometry() == 0, "viewport image rebuilt");
    CHECK(viewport->getHUDKit() != NULL, "viewport image owns HUD kit");
    CHECK(viewport->getTextureNode() != NULL, "viewport image owns texture node");
    CHECK(viewport->getImageFaceSet() != NULL, "viewport image owns face set");
    CHECK(texture_matches(viewport->getTextureNode(), 4, 2, 4),
	  "viewport texture dimensions match source");
    CHECK(viewport->realizedDataRevision.getValue() == source->dataRevision.getValue(),
	  "viewport data revision recorded");
    CHECK(viewport->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "viewport dirty revision recorded");

    unsigned char pixel[4] = {200, 30, 20, 255};
    CHECK(imgstream_write_rect(stream, 1, 0, 1, 1, pixel, 4) == 0,
	  "viewport source dirty write accepted");
    CHECK(viewport->syncFromSource() == 0, "viewport sync rebuilt after dirty source");
    CHECK(viewport->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "viewport dirty revision refreshed");
    CHECK(texture_matches(viewport->getTextureNode(), 4, 2, 4),
	  "viewport refreshed texture dimensions match source");

    viewport->unref();
    source->unref();
    imgstream_destroy(stream);
    return 0;
}

static int
test_image_plane(void)
{
    imgstream_t *stream = make_rgba_stream();
    CHECK(stream != NULL, "created RGBA stream for image plane");

    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();
    CHECK(source->setStream(stream) == 0, "attached plane source stream");

    SoBRLImagePlane *plane = new SoBRLImagePlane;
    plane->ref();
    plane->imageSource.setValue(source);
    plane->sourcePath = "/render/image";
    plane->sourceId = 42;
    plane->sizeMode = SoBRLImagePlane::EXPLICIT_SIZE;
    plane->width = 2.0f;
    plane->height = 1.0f;
    plane->fit = SoBRLImagePlane::STRETCH;
    plane->selectable = FALSE;
    plane->depthTest = TRUE;
    plane->depthWrite = FALSE;
    plane->doubleSided = TRUE;

    CHECK(plane->rebuildGeometry() == 0, "image plane rebuilt");
    CHECK(plane->getNumChildren() == 1, "image plane owns one textured quad");
    CHECK(plane->getTextureNode() != NULL, "image plane owns texture node");
    CHECK(plane->getImageFaceSet() != NULL, "image plane owns face set");
    CHECK(texture_matches(plane->getTextureNode(), 4, 2, 4),
	  "image plane texture dimensions match source");
    CHECK(plane->realizedDataRevision.getValue() == source->dataRevision.getValue(),
	  "image plane data revision recorded");
    CHECK(plane->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "image plane dirty revision recorded");

    unsigned char pixel[4] = {0, 180, 60, 255};
    CHECK(imgstream_write_rect(stream, 2, 1, 1, 1, pixel, 4) == 0,
	  "plane source dirty write accepted");
    CHECK(plane->syncFromSource() == 0, "image plane sync rebuilt after dirty source");
    CHECK(plane->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	  "image plane dirty revision refreshed");

    plane->unref();
    source->unref();
    imgstream_destroy(stream);
    return 0;
}

int
main(int ac, char **av)
{
    bu_setprogname(av[0]);
    (void)ac;
    (void)av;

    bobol_init(NULL);

    if (test_viewport_image())
	return 1;
    if (test_image_plane())
	return 1;

    return 0;
}
