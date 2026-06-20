/*              I M A G E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/image_source.h"

#include <limits.h>
#include <stdint.h>

SO_NODE_SOURCE(SoBRLImageSource);

static uint32_t
image_source_clamp_u32(uint64_t value)
{
    return value > UINT32_MAX ? UINT32_MAX : (uint32_t)value;
}

static uint32_t
next_revision(uint32_t value)
{
    return value == UINT32_MAX ? UINT32_MAX : value + 1;
}

static SoBRLImageSource::PixelFormat
pixel_format_from_stream(enum imgstream_pixel_format format)
{
    switch (format) {
	case IMGSTREAM_PIXEL_RGB8:
	    return SoBRLImageSource::PIXEL_RGB8;
	case IMGSTREAM_PIXEL_RGBA8:
	    return SoBRLImageSource::PIXEL_RGBA8;
	default:
	    return SoBRLImageSource::PIXEL_UNKNOWN;
    }
}

SoBRLImageSource::SoBRLImageSource(void) :
    stream(NULL),
    subscriberId(-1),
    streamOwned(FALSE)
{
    SO_NODE_CONSTRUCTOR(SoBRLImageSource);

    SO_NODE_DEFINE_ENUM_VALUE(SourceKind, SOURCE_EMPTY);
    SO_NODE_DEFINE_ENUM_VALUE(SourceKind, SOURCE_STATIC_IMAGE);
    SO_NODE_DEFINE_ENUM_VALUE(SourceKind, SOURCE_IMAGE_STREAM);
    SO_NODE_DEFINE_ENUM_VALUE(Status, STATUS_EMPTY);
    SO_NODE_DEFINE_ENUM_VALUE(Status, STATUS_READY);
    SO_NODE_DEFINE_ENUM_VALUE(Status, STATUS_STREAMING);
    SO_NODE_DEFINE_ENUM_VALUE(Status, STATUS_FAILED);
    SO_NODE_DEFINE_ENUM_VALUE(PixelFormat, PIXEL_UNKNOWN);
    SO_NODE_DEFINE_ENUM_VALUE(PixelFormat, PIXEL_RGB8);
    SO_NODE_DEFINE_ENUM_VALUE(PixelFormat, PIXEL_RGBA8);

    SO_NODE_ADD_FIELD(imageId, (""));
    SO_NODE_ADD_FIELD(sourceUri, (""));
    SO_NODE_ADD_FIELD(sourceKind, (SOURCE_EMPTY));
    SO_NODE_SET_SF_ENUM_TYPE(sourceKind, SourceKind);
    SO_NODE_ADD_FIELD(status, (STATUS_EMPTY));
    SO_NODE_SET_SF_ENUM_TYPE(status, Status);
    SO_NODE_ADD_FIELD(diagnostic, (""));
    SO_NODE_ADD_FIELD(pixelWidth, (0));
    SO_NODE_ADD_FIELD(pixelHeight, (0));
    SO_NODE_ADD_FIELD(pixelFormat, (PIXEL_UNKNOWN));
    SO_NODE_SET_SF_ENUM_TYPE(pixelFormat, PixelFormat);
    SO_NODE_ADD_FIELD(colorSpace, (""));
    SO_NODE_ADD_FIELD(alphaMode, ("none"));
    SO_NODE_ADD_FIELD(sourceRevision, (0));
    SO_NODE_ADD_FIELD(dataRevision, (0));
    SO_NODE_ADD_FIELD(dirtyRevision, (0));
    SO_NODE_ADD_FIELD(dirtyX, (0));
    SO_NODE_ADD_FIELD(dirtyY, (0));
    SO_NODE_ADD_FIELD(dirtyWidth, (0));
    SO_NODE_ADD_FIELD(dirtyHeight, (0));
    SO_NODE_ADD_FIELD(streamConnected, (FALSE));
    SO_NODE_ADD_FIELD(producerActive, (FALSE));
}

SoBRLImageSource::~SoBRLImageSource(void)
{
    this->releaseStream();
}

void
SoBRLImageSource::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLImageSource, SoNode, "Node");
}

void
SoBRLImageSource::resetStreamFields(Status statusValue, const char *diagnosticValue)
{
    this->sourceKind = SOURCE_EMPTY;
    this->status = statusValue;
    this->diagnostic = diagnosticValue ? diagnosticValue : "";
    this->pixelWidth = 0;
    this->pixelHeight = 0;
    this->pixelFormat = PIXEL_UNKNOWN;
    this->colorSpace = "";
    this->alphaMode = "none";
    this->dataRevision = 0;
    this->dirtyRevision = 0;
    this->dirtyX = 0;
    this->dirtyY = 0;
    this->dirtyWidth = 0;
    this->dirtyHeight = 0;
    this->streamConnected = FALSE;
    this->producerActive = FALSE;
}

void
SoBRLImageSource::releaseStream(void)
{
    if (this->stream && this->subscriberId >= 0)
	imgstream_unsubscribe(this->stream, this->subscriberId);

    if (this->stream && this->streamOwned)
	imgstream_destroy(this->stream);

    this->stream = NULL;
    this->subscriberId = -1;
    this->streamOwned = FALSE;
}

void
SoBRLImageSource::clearSource(void)
{
    this->releaseStream();
    this->sourceUri = "";
    this->resetStreamFields(STATUS_EMPTY, "");
    this->sourceRevision = next_revision(this->sourceRevision.getValue());
}

imgstream_t *
SoBRLImageSource::getStream(void) const
{
    return this->stream;
}

SbBool
SoBRLImageSource::ownsStream(void) const
{
    return this->streamOwned;
}

int
SoBRLImageSource::setStream(imgstream_t *newStream)
{
    if (!newStream) {
	this->clearSource();
	this->status = STATUS_FAILED;
	this->diagnostic = "image stream is null";
	return -1;
    }

    return this->attachStream(newStream, FALSE, SOURCE_IMAGE_STREAM);
}

int
SoBRLImageSource::setImage(const icv_image_t *image)
{
    if (!image) {
	this->clearSource();
	this->status = STATUS_FAILED;
	this->diagnostic = "image is null";
	return -1;
    }

    imgstream_t *newStream = imgstream_create_from_icv(image);
    if (!newStream) {
	this->clearSource();
	this->status = STATUS_FAILED;
	this->diagnostic = "failed to create stream from image";
	return -1;
    }

    return this->attachStream(newStream, TRUE, SOURCE_STATIC_IMAGE);
}

int
SoBRLImageSource::attachStream(imgstream_t *newStream, SbBool owned, SourceKind kind)
{
    this->releaseStream();

    this->stream = newStream;
    this->streamOwned = owned;
    this->subscriberId = imgstream_subscribe(this->stream, SoBRLImageSource::dirtyCB, this);
    if (this->subscriberId < 0) {
	if (owned)
	    imgstream_destroy(newStream);
	this->stream = NULL;
	this->streamOwned = FALSE;
	this->resetStreamFields(STATUS_FAILED, "failed to subscribe to image stream");
	this->sourceRevision = next_revision(this->sourceRevision.getValue());
	return -1;
    }

    this->sourceKind = kind;
    this->streamConnected = TRUE;
    this->sourceUri = owned ? "icv:memory" : "";
    this->sourceRevision = next_revision(this->sourceRevision.getValue());
    return this->refreshFromStream();
}

int
SoBRLImageSource::refreshFromStream(void)
{
    if (!this->stream) {
	this->resetStreamFields(STATUS_EMPTY, "");
	return -1;
    }

    struct imgstream_info info;
    if (imgstream_get_info(this->stream, &info) != 0) {
	this->status = STATUS_FAILED;
	this->diagnostic = "failed to query image stream";
	return -1;
    }

    this->pixelWidth = image_source_clamp_u32(info.width);
    this->pixelHeight = image_source_clamp_u32(info.height);
    this->pixelFormat = pixel_format_from_stream(info.format);
    this->colorSpace = "srgb";
    this->alphaMode = info.format == IMGSTREAM_PIXEL_RGBA8 ? "straight" : "none";
    this->dataRevision = image_source_clamp_u32(info.generation);
    this->producerActive = info.producer_active ? TRUE : FALSE;
    this->status = info.producer_active ? STATUS_STREAMING : STATUS_READY;
    this->diagnostic = "";
    this->streamConnected = TRUE;

    if (info.dirty) {
	this->dirtyRevision = image_source_clamp_u32(info.generation);
	this->dirtyX = image_source_clamp_u32(info.dirty_rect.x);
	this->dirtyY = image_source_clamp_u32(info.dirty_rect.y);
	this->dirtyWidth = image_source_clamp_u32(info.dirty_rect.width);
	this->dirtyHeight = image_source_clamp_u32(info.dirty_rect.height);
    }

    return 0;
}

void
SoBRLImageSource::dirtyCB(void *ctx, const struct imgstream_rect *rect, uint64_t generation)
{
    SoBRLImageSource *source = static_cast<SoBRLImageSource *>(ctx);
    if (!source)
	return;

    source->dataRevision = image_source_clamp_u32(generation);
    source->dirtyRevision = image_source_clamp_u32(generation);
    if (rect) {
	source->dirtyX = image_source_clamp_u32(rect->x);
	source->dirtyY = image_source_clamp_u32(rect->y);
	source->dirtyWidth = image_source_clamp_u32(rect->width);
	source->dirtyHeight = image_source_clamp_u32(rect->height);
    }

    source->refreshFromStream();
}
