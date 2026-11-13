/*              I M A G E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BImageSource.h"
#include "image_display_util.h"
#include "image_source_private.h"
#include "identity_counter_private.h"
#include "scalar_publication_private.h"

#include <limits.h>
#include <stdint.h>

#include <exception>
#include <memory>
#include <stdexcept>

SO_NODE_SOURCE(SoBRLImageSource);

static const char image_source_alpha_none[] = "none";
static const char image_source_alpha_straight[] = "straight";
static const char image_source_color_space[] = "srgb";
static const char image_source_memory_uri[] = "icv:memory";
static const char image_source_null_image[] = "image is null";
static const char image_source_null_stream[] = "image stream is null";
static const char image_source_stream_create_failure[] =
    "failed to create stream from image";
static const char image_source_stream_query_failure[] =
    "failed to query image stream";
static const char image_source_stream_subscribe_failure[] =
    "failed to subscribe to image stream";

static uint32_t
image_source_clamp_u32(uint64_t value)
{
    return value > UINT32_MAX ? UINT32_MAX : (uint32_t)value;
}

static uint32_t
next_revision(uint32_t value)
{
    return bobol_identity_successor_or_terminate(value);
}

static void
atomic_store_max(std::atomic<uint64_t> &target, uint64_t value)
{
    uint64_t current = target.load(std::memory_order_acquire);
    while (current < value &&
	    !target.compare_exchange_weak(current, value,
		std::memory_order_release, std::memory_order_acquire))
	;
}

struct BObolImageSourceSubscription {
    BObolImageSourceSubscription(imgstream_t *sourceStream, SbBool takeOwnership,
	SoBRLImageSource::SourceKind sourceKind) :
	stream(sourceStream), subscriberId(-1), owned(takeOwnership), kind(sourceKind)
    {
    }

    ~BObolImageSourceSubscription()
    {
	target.store(nullptr);
	if (stream && subscriberId >= 0)
	    imgstream_unsubscribe(stream, subscriberId);
	if (stream && owned)
	    imgstream_destroy(stream);
    }

    imgstream_t *stream;
    imgstream_subscriber_id subscriberId;
    SbBool owned;
    SoBRLImageSource::SourceKind kind;
    std::atomic<uint64_t> pendingGeneration{0};
    std::atomic<SoBRLImageSource *> target{nullptr};
};

struct BObolImageSourcePublication::Impl {
    explicit Impl(SoBRLImageSource &target) :
	node(target), owner(new SoBRLImageSource),
	candidate(*static_cast<SoBRLImageSource *>(owner.get())),
	nextRealized(target.realizedGeneration.load(std::memory_order_acquire))
    {
	copy_publication_scalar_fields(candidate, node);
    }

    SoBRLImageSource &node;
    SbModernUtils::SoNodeRef owner;
    SoBRLImageSource &candidate;
    std::unique_ptr<BObolImageSourceSubscription> nextSubscription;
    std::unique_ptr<PreparedScalarFields> fields;
    struct imgstream_info refreshInfo = {};
    uint64_t nextRealized;
    bool replace = false;
    bool hasRefreshInfo = false;
    bool committed = false;
};

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

static void
image_source_reset_fields(SoBRLImageSource &source,
	SoBRLImageSource::Status status, const char *diagnostic)
{
    source.sourceUri = "";
    source.sourceKind = SoBRLImageSource::SOURCE_EMPTY;
    source.status = status;
    source.diagnostic = diagnostic ? diagnostic : "";
    source.pixelWidth = 0;
    source.pixelHeight = 0;
    source.pixelFormat = SoBRLImageSource::PIXEL_UNKNOWN;
    source.colorSpace = "";
    source.alphaMode = image_source_alpha_none;
    source.dataRevision = 0;
    source.dirtyRevision = 0;
    source.dirtyX = 0;
    source.dirtyY = 0;
    source.dirtyWidth = 0;
    source.dirtyHeight = 0;
    source.streamConnected = FALSE;
    source.producerActive = FALSE;
}

static bool
image_source_has_reset_fields(const SoBRLImageSource &source,
	SoBRLImageSource::Status status, const char *diagnostic)
{
    const char *message = diagnostic ? diagnostic : "";
    return source.sourceUri.getValue().getLength() == 0 &&
	source.sourceKind.getValue() == SoBRLImageSource::SOURCE_EMPTY &&
	source.status.getValue() == status &&
	source.diagnostic.getValue() == message &&
	!source.pixelWidth.getValue() && !source.pixelHeight.getValue() &&
	source.pixelFormat.getValue() == SoBRLImageSource::PIXEL_UNKNOWN &&
	source.colorSpace.getValue().getLength() == 0 &&
	source.alphaMode.getValue() == image_source_alpha_none &&
	!source.dataRevision.getValue() && !source.dirtyRevision.getValue() &&
	!source.dirtyX.getValue() && !source.dirtyY.getValue() &&
	!source.dirtyWidth.getValue() && !source.dirtyHeight.getValue() &&
	!source.streamConnected.getValue() && !source.producerActive.getValue();
}

static int
image_source_publish_failure(SoBRLImageSource &source, const char *diagnostic)
{
    if (!source.getStream() &&
	image_source_has_reset_fields(source, SoBRLImageSource::STATUS_FAILED,
	    diagnostic))
	return -1;

    BObolImageSourcePublication publication(source);
    image_source_reset_fields(publication.next(),
	SoBRLImageSource::STATUS_FAILED, diagnostic);
    publication.next().sourceRevision =
	next_revision(source.sourceRevision.getValue());
    publication.replaceSubscription(nullptr, 0);
    publication.publish();
    return -1;
}

static void
image_source_apply_info(SoBRLImageSource &source,
	const BObolImageSourceSubscription &subscription,
	const struct imgstream_info &info)
{
    source.sourceUri = subscription.owned ? image_source_memory_uri : "";
    source.sourceKind = subscription.kind;
    source.pixelWidth = image_source_clamp_u32(info.width);
    source.pixelHeight = image_source_clamp_u32(info.height);
    source.pixelFormat = pixel_format_from_stream(info.format);
    source.colorSpace = image_source_color_space;
    source.alphaMode = info.format == IMGSTREAM_PIXEL_RGBA8 ?
	image_source_alpha_straight : image_source_alpha_none;
    source.dataRevision = image_source_clamp_u32(info.generation);
    source.producerActive = info.producer_active ? TRUE : FALSE;
    source.status = info.producer_active ? SoBRLImageSource::STATUS_STREAMING :
	SoBRLImageSource::STATUS_READY;
    source.diagnostic = "";
    source.streamConnected = TRUE;
    if (info.dirty) {
	source.dirtyRevision = image_source_clamp_u32(info.generation);
	source.dirtyX = image_source_clamp_u32(info.dirty_rect.x);
	source.dirtyY = image_source_clamp_u32(info.dirty_rect.y);
	source.dirtyWidth = image_source_clamp_u32(info.dirty_rect.width);
	source.dirtyHeight = image_source_clamp_u32(info.dirty_rect.height);
    }
}

static bool
image_source_matches_info(const SoBRLImageSource &source,
	const BObolImageSourceSubscription &subscription,
	const struct imgstream_info &info)
{
    if (source.sourceUri.getValue() !=
	    (subscription.owned ? image_source_memory_uri : "") ||
	source.sourceKind.getValue() != subscription.kind ||
	source.pixelWidth.getValue() != image_source_clamp_u32(info.width) ||
	source.pixelHeight.getValue() != image_source_clamp_u32(info.height) ||
	source.pixelFormat.getValue() != pixel_format_from_stream(info.format) ||
	source.colorSpace.getValue() != image_source_color_space ||
	source.alphaMode.getValue() !=
	    (info.format == IMGSTREAM_PIXEL_RGBA8 ?
		image_source_alpha_straight : image_source_alpha_none) ||
	source.dataRevision.getValue() != image_source_clamp_u32(info.generation) ||
	source.producerActive.getValue() != (info.producer_active ? TRUE : FALSE) ||
	source.status.getValue() != (info.producer_active ?
	    SoBRLImageSource::STATUS_STREAMING : SoBRLImageSource::STATUS_READY) ||
	source.diagnostic.getValue().getLength() != 0 ||
	!source.streamConnected.getValue())
	return false;
    return !info.dirty ||
	(source.dirtyRevision.getValue() ==
	    image_source_clamp_u32(info.generation) &&
	 source.dirtyX.getValue() == image_source_clamp_u32(info.dirty_rect.x) &&
	 source.dirtyY.getValue() == image_source_clamp_u32(info.dirty_rect.y) &&
	 source.dirtyWidth.getValue() ==
	    image_source_clamp_u32(info.dirty_rect.width) &&
	 source.dirtyHeight.getValue() ==
	    image_source_clamp_u32(info.dirty_rect.height));
}

BObolImageSourcePublication::BObolImageSourcePublication(
    SoBRLImageSource &target) :
    impl(new Impl(target))
{
}

BObolImageSourcePublication::~BObolImageSourcePublication() = default;

int
BObolImageSourcePublication::adoptStream(SoBRLImageSource &target,
    imgstream_t *stream)
{
    return target.attachStream(stream, TRUE,
	SoBRLImageSource::SOURCE_IMAGE_STREAM);
}

BObolImageSourceRefresh
BObolImageSourcePublication::queryRefresh(SoBRLImageSource &target,
    struct imgstream_info &info)
{
    if (!target.subscription ||
	imgstream_get_info(target.subscription->stream, &info) != 0)
	return BObolImageSourceRefresh::Failed;

    atomic_store_max(target.pendingGeneration, info.generation);
    const bool current =
	image_source_matches_info(target, *target.subscription, info) &&
	target.realizedGeneration.load(std::memory_order_acquire) ==
	    info.generation;
    return current ? BObolImageSourceRefresh::Current :
	BObolImageSourceRefresh::Required;
}

BObolImageSourceRefresh
BObolImageSourcePublication::refreshRequired(SoBRLImageSource &target)
{
    struct imgstream_info info;
    return queryRefresh(target, info);
}

SoBRLImageSource &
BObolImageSourcePublication::next()
{
    if (this->impl->fields || this->impl->committed)
	throw std::logic_error("image source publication already prepared");
    return this->impl->candidate;
}

const SoBRLImageSource &
BObolImageSourcePublication::successor() const
{
    return this->impl->candidate;
}

void
BObolImageSourcePublication::replaceSubscription(
    std::unique_ptr<BObolImageSourceSubscription> replacement,
    uint64_t realizedGeneration)
{
    if (this->impl->fields || this->impl->committed)
	throw std::logic_error("image source publication already prepared");
    this->impl->nextSubscription = std::move(replacement);
    this->impl->nextRealized = realizedGeneration;
    this->impl->replace = true;
}

void
BObolImageSourcePublication::setRealizedGeneration(uint64_t generation)
{
    if (this->impl->fields || this->impl->committed)
	throw std::logic_error("image source publication already prepared");
    this->impl->nextRealized = generation;
}

BObolImageSourceRefresh
BObolImageSourcePublication::prepareRefresh()
{
    if (this->impl->fields || this->impl->committed)
	throw std::logic_error("image source publication already prepared");
    const BObolImageSourceRefresh state = queryRefresh(this->impl->node,
	this->impl->refreshInfo);
    if (state == BObolImageSourceRefresh::Failed)
	return state;
    this->impl->hasRefreshInfo = true;
    if (state == BObolImageSourceRefresh::Current)
	return state;

    image_source_apply_info(this->impl->candidate,
	*this->impl->node.subscription, this->impl->refreshInfo);
    this->impl->nextRealized = this->impl->refreshInfo.generation;
    this->prepare();
    return BObolImageSourceRefresh::Required;
}

int
BObolImageSourcePublication::loadPreparedPayload(
    struct bobol_image_payload &payload) const
{
    if (!this->impl->hasRefreshInfo || !this->impl->node.subscription)
	return -1;
    return bobol_image_payload_load_info(
	this->impl->node.subscription->stream, this->impl->refreshInfo,
	image_source_clamp_u32(this->impl->candidate.dirtyRevision.getValue()),
	&payload);
}

void
BObolImageSourcePublication::prepare()
{
    if (this->impl->fields || this->impl->committed)
	throw std::logic_error("image source publication already prepared");
    this->impl->fields = std::make_unique<PreparedScalarFields>(
	this->impl->node, this->impl->candidate);
}

void
BObolImageSourcePublication::commit()
{
    if (!this->impl->fields || this->impl->committed)
	throw std::logic_error("invalid image source publication commit");

    BObolImageSourceSubscription *previous = nullptr;
    if (this->impl->replace) {
	previous = this->impl->node.subscription;
	if (previous)
	    previous->target.store(nullptr);
	this->impl->node.subscription =
	    this->impl->nextSubscription.release();
	this->impl->node.subscriberId = this->impl->node.subscription ?
	    this->impl->node.subscription->subscriberId : -1;
	this->impl->node.streamOwned = this->impl->node.subscription ?
	    this->impl->node.subscription->owned : FALSE;
	this->impl->node.pendingGeneration.store(this->impl->nextRealized,
	    std::memory_order_release);
	if (this->impl->node.subscription) {
	    this->impl->node.subscription->target.store(&this->impl->node);
	    atomic_store_max(this->impl->node.pendingGeneration,
		this->impl->node.subscription->pendingGeneration.load());
	}
    }
    this->impl->node.realizedGeneration.store(this->impl->nextRealized,
	std::memory_order_release);
    this->impl->fields->commit();
    this->impl->committed = true;
    delete previous;
}

void
BObolImageSourcePublication::restore()
{
    if (this->impl->fields)
	this->impl->fields->restore();
}

void
BObolImageSourcePublication::notify(std::exception_ptr &failure)
{
    if (!this->impl->committed)
	return;
    this->restore();
    if (this->impl->fields->changed())
	this->impl->fields->notify(failure);
}

void
BObolImageSourcePublication::notify()
{
    std::exception_ptr failure;
    this->notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}

void
BObolImageSourcePublication::publish()
{
    this->prepare();
    this->commit();
    this->notify();
}

SoBRLImageSource::SoBRLImageSource(void) :
    subscription(nullptr),
    subscriberId(-1),
    pendingGeneration(0),
    realizedGeneration(0),
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
    SO_NODE_ADD_FIELD(alphaMode, (image_source_alpha_none));
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
SoBRLImageSource::releaseStream(void)
{
    delete this->subscription;
    this->subscription = nullptr;
    this->subscriberId = -1;
    this->streamOwned = FALSE;
    this->pendingGeneration.store(0, std::memory_order_release);
    this->realizedGeneration.store(0, std::memory_order_release);
}

void
SoBRLImageSource::clearSource(void)
{
    if (!this->subscription &&
	image_source_has_reset_fields(*this, STATUS_EMPTY, ""))
	return;
    BObolImageSourcePublication publication(*this);
    image_source_reset_fields(publication.next(), STATUS_EMPTY, "");
    publication.next().sourceRevision =
	next_revision(this->sourceRevision.getValue());
    publication.replaceSubscription(nullptr, 0);
    publication.publish();
}

imgstream_t *
SoBRLImageSource::getStream(void) const
{
    return this->subscription ? this->subscription->stream : nullptr;
}

SbBool
SoBRLImageSource::ownsStream(void) const
{
    return this->streamOwned;
}

SbBool
SoBRLImageSource::hasPendingStreamUpdate(void) const
{
    uint64_t pending = this->pendingGeneration.load(std::memory_order_acquire);
    uint64_t realized = this->realizedGeneration.load(std::memory_order_acquire);
    return pending > realized ? TRUE : FALSE;
}

int
SoBRLImageSource::setStream(imgstream_t *newStream)
{
    if (!newStream)
	return image_source_publish_failure(*this, image_source_null_stream);

    if (newStream == this->getStream())
	return 0;

    return this->attachStream(newStream, FALSE, SOURCE_IMAGE_STREAM);
}

int
SoBRLImageSource::setImage(const icv_image_t *image)
{
    if (!image)
	return image_source_publish_failure(*this, image_source_null_image);

    imgstream_t *newStream = imgstream_create_from_icv(image);
    if (!newStream)
	return image_source_publish_failure(*this,
	    image_source_stream_create_failure);

    return this->attachStream(newStream, TRUE, SOURCE_STATIC_IMAGE);
}

int
SoBRLImageSource::attachStream(imgstream_t *newStream, SbBool owned, SourceKind kind)
{
    using StreamOwner = std::unique_ptr<imgstream_t, void (*)(imgstream_t *)>;
    StreamOwner ownedStream(owned ? newStream : nullptr, imgstream_destroy);
    BObolImageSourcePublication publication(*this);
    auto nextSubscription = std::make_unique<BObolImageSourceSubscription>(
	newStream, owned, kind);
    if (owned)
	ownedStream.release();

    nextSubscription->subscriberId = imgstream_subscribe(newStream,
	SoBRLImageSource::dirtyCB, nextSubscription.get());
    if (nextSubscription->subscriberId < 0) {
	image_source_reset_fields(publication.next(), STATUS_FAILED,
	    image_source_stream_subscribe_failure);
	publication.next().sourceRevision =
	    next_revision(this->sourceRevision.getValue());
	publication.replaceSubscription(nullptr, 0);
	publication.publish();
	return -1;
    }

    image_source_reset_fields(publication.next(), STATUS_EMPTY, "");
    publication.next().sourceKind = kind;
    publication.next().sourceUri = owned ? image_source_memory_uri : "";
    publication.next().streamConnected = TRUE;
    publication.next().sourceRevision =
	next_revision(this->sourceRevision.getValue());

    struct imgstream_info info;
    const bool valid = imgstream_get_info(newStream, &info) == 0;
    uint64_t realized = 0;
    if (valid) {
	image_source_apply_info(publication.next(), *nextSubscription, info);
	nextSubscription->pendingGeneration.store(info.generation);
	realized = info.generation;
    } else {
	publication.next().status = STATUS_FAILED;
	publication.next().diagnostic = image_source_stream_query_failure;
    }
    publication.replaceSubscription(std::move(nextSubscription), realized);
    publication.publish();
    return valid ? 0 : -1;
}

int
SoBRLImageSource::refreshFromStream(void)
{
    if (!this->subscription) {
	if (!image_source_has_reset_fields(*this, STATUS_EMPTY, "")) {
	    BObolImageSourcePublication publication(*this);
	    image_source_reset_fields(publication.next(), STATUS_EMPTY, "");
	    publication.publish();
	}
	return -1;
    }

    struct imgstream_info info;
    if (imgstream_get_info(this->subscription->stream, &info) != 0) {
	if (this->status.getValue() != STATUS_FAILED ||
	    this->diagnostic.getValue() != image_source_stream_query_failure) {
	    BObolImageSourcePublication publication(*this);
	    publication.next().sourceUri =
		this->subscription->owned ? image_source_memory_uri : "";
	    publication.next().sourceKind = this->subscription->kind;
	    publication.next().status = STATUS_FAILED;
	    publication.next().diagnostic = image_source_stream_query_failure;
	    publication.next().streamConnected = TRUE;
	    publication.publish();
	}
	return -1;
    }

    atomic_store_max(this->pendingGeneration, info.generation);
    if (image_source_matches_info(*this, *this->subscription, info) &&
	this->realizedGeneration.load(std::memory_order_acquire) == info.generation)
	return 0;

    BObolImageSourcePublication publication(*this);
    image_source_apply_info(publication.next(), *this->subscription, info);
    publication.setRealizedGeneration(info.generation);
    publication.publish();

    return 0;
}

void
SoBRLImageSource::dirtyCB(void *ctx, const struct imgstream_rect *UNUSED(rect), uint64_t generation)
{
    auto *subscription = static_cast<BObolImageSourceSubscription *>(ctx);
    if (!subscription)
	return;

    atomic_store_max(subscription->pendingGeneration, generation);
    SoBRLImageSource *source = subscription->target.load();
    if (source)
	atomic_store_max(source->pendingGeneration, generation);
}
