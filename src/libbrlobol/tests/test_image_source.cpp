/*              T E S T _ I M A G E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "brlobol.h"

#include "bu/log.h"
#include "bu/str.h"
#include "icv.h"
#include "imgstream/stream.h"

#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

static int
test_external_stream_source(void)
{
    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();

    CHECK(source->sourceKind.getValue() == SoBRLImageSource::SOURCE_EMPTY,
	  "default source kind is empty");
    CHECK(source->status.getValue() == SoBRLImageSource::STATUS_EMPTY,
	  "default status is empty");
    CHECK(source->pixelWidth.getValue() == 0 && source->pixelHeight.getValue() == 0,
	  "default dimensions are empty");
    CHECK(source->getStream() == NULL, "default stream is null");

    imgstream_t *stream = imgstream_create(4, 3, IMGSTREAM_PIXEL_RGB8);
    CHECK(stream != NULL, "created external stream");
    CHECK(source->setStream(stream) == 0, "external stream attached");
    CHECK(source->getStream() == stream, "external stream stored");
    CHECK(source->ownsStream() == FALSE, "external stream remains borrowed");
    CHECK(source->streamConnected.getValue() == TRUE, "stream marked connected");
    CHECK(source->sourceKind.getValue() == SoBRLImageSource::SOURCE_IMAGE_STREAM,
	  "external source kind recorded");
    CHECK(source->status.getValue() == SoBRLImageSource::STATUS_READY,
	  "external stream is initially ready");
    CHECK(source->pixelWidth.getValue() == 4 && source->pixelHeight.getValue() == 3,
	  "external stream dimensions recorded");
    CHECK(source->pixelFormat.getValue() == SoBRLImageSource::PIXEL_RGB8,
	  "external stream pixel format recorded");

    unsigned char pixels[2 * 2 * 3];
    memset(pixels, 127, sizeof(pixels));
    CHECK(imgstream_write_rect(stream, 1, 1, 2, 2, pixels, 2 * 3) == 0,
	  "external stream rect write accepted");
    CHECK(source->hasPendingStreamUpdate() == TRUE,
	  "rect write marked stream update pending");
    CHECK(source->dataRevision.getValue() == 0,
	  "rect write did not mutate data revision before refresh");
    CHECK(source->dirtyRevision.getValue() == 0,
	  "rect write did not mutate dirty revision before refresh");
    CHECK(source->refreshFromStream() == 0, "pending rect write refresh accepted");
    CHECK(source->dataRevision.getValue() == 1, "refresh updated data revision");
    CHECK(source->dirtyRevision.getValue() == 1, "refresh updated dirty revision");
    CHECK(source->dirtyX.getValue() == 1 && source->dirtyY.getValue() == 1,
	  "rect write dirty origin recorded");
    CHECK(source->dirtyWidth.getValue() == 2 && source->dirtyHeight.getValue() == 2,
	  "rect write dirty size recorded");
    CHECK(source->hasPendingStreamUpdate() == FALSE,
	  "refresh consumed pending stream update");

    CHECK(imgstream_producer_begin(stream) == 0, "producer begin accepted");
    CHECK(source->refreshFromStream() == 0, "active producer refresh accepted");
    CHECK(source->producerActive.getValue() == TRUE, "producer active recorded");
    CHECK(source->status.getValue() == SoBRLImageSource::STATUS_STREAMING,
	  "active producer status recorded");
    CHECK(imgstream_producer_end(stream) == 0, "producer end accepted");
    CHECK(source->refreshFromStream() == 0, "inactive producer refresh accepted");
    CHECK(source->producerActive.getValue() == FALSE, "inactive producer recorded");
    CHECK(source->status.getValue() == SoBRLImageSource::STATUS_READY,
	  "inactive producer status recorded");

    source->clearSource();
    CHECK(source->getStream() == NULL, "clear released external stream");
    CHECK(source->streamConnected.getValue() == FALSE, "clear disconnected source");
    CHECK(source->dirtyRevision.getValue() == 0, "clear reset dirty revision");
    CHECK(imgstream_write_rect(stream, 0, 0, 1, 1, pixels, 1 * 3) == 0,
	  "write after clear still accepted by external owner");
    CHECK(source->dirtyRevision.getValue() == 0,
	  "cleared source no longer receives external stream notifications");
    CHECK(source->hasPendingStreamUpdate() == FALSE,
	  "cleared source ignores external stream dirty markers");

    imgstream_destroy(stream);
    source->unref();
    return 0;
}

static int
test_static_image_source(void)
{
    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();

    icv_image_t *image = icv_create(2, 2, ICV_COLOR_SPACE_RGB);
    CHECK(image != NULL, "created ICV source image");
    for (size_t i = 0; i < image->width * image->height * (size_t)image->channels; i++)
	image->data[i] = (double)(i % 3) / 2.0;

    CHECK(source->setImage(image) == 0, "static ICV image attached");
    icv_destroy(image);

    CHECK(source->getStream() != NULL, "static image created stream");
    CHECK(source->ownsStream() == TRUE, "static image stream is owned");
    CHECK(source->sourceKind.getValue() == SoBRLImageSource::SOURCE_STATIC_IMAGE,
	  "static source kind recorded");
    CHECK(source->status.getValue() == SoBRLImageSource::STATUS_READY,
	  "static source is ready");
    CHECK(source->pixelWidth.getValue() == 2 && source->pixelHeight.getValue() == 2,
	  "static dimensions recorded");
    CHECK(source->pixelFormat.getValue() == SoBRLImageSource::PIXEL_RGB8,
	  "static image pixel format recorded");
    CHECK(source->dataRevision.getValue() == 1, "static image copied at generation one");
    CHECK(source->dirtyRevision.getValue() == 1, "static image initial dirty generation recorded");
    CHECK(BU_STR_EQUAL(source->colorSpace.getValue().getString(), "srgb"),
	  "static image color space recorded");
    CHECK(BU_STR_EQUAL(source->alphaMode.getValue().getString(), "none"),
	  "static image alpha mode recorded");

    source->clearSource();
    CHECK(source->getStream() == NULL, "clear destroyed owned static stream");
    CHECK(source->ownsStream() == FALSE, "clear reset owned flag");

    source->unref();
    return 0;
}

int
main(int ac, char **av)
{
    bu_setprogname(av[0]);
    (void)ac;
    (void)av;

    brlobol_init(NULL);

    if (test_external_stream_source())
	return 1;
    if (test_static_image_source())
	return 1;

    return 0;
}
