/*          I M A G E _ D I S P L A Y _ U T I L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BViewController.h"
#include "BObol/BViewportImage.h"
#include "image_display_util.h"
#include "scalar_publication_private.h"

#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPickStyle.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoTextureCoordinate2.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>

#include <algorithm>
#include <array>
#include <exception>
#include <limits.h>
#include <stdexcept>

struct BObolPreparedViewportImage::Impl {
    explicit Impl(SoBRLViewportImage &viewport) :
	target(viewport), candidateOwner(new SoBRLViewportImage)
    {
    }

    SoBRLViewportImage &target;
    SbModernUtils::SoNodeRef candidateOwner;
    std::unique_ptr<PreparedScalarFields> fields;
    std::unique_ptr<SoChildList::Replacement> children;
    SoTexture2 *texture = NULL;
    SoFaceSet *face = NULL;
    bool ready = false;
    bool committed = false;
};

BObolPreparedViewportImage::BObolPreparedViewportImage(
    SoBRLViewportImage &target) :
    impl(new Impl(target))
{
    auto *candidate =
	static_cast<SoBRLViewportImage *>(this->impl->candidateOwner.get());
    copy_publication_scalar_fields(*candidate, target);
}

BObolPreparedViewportImage::BObolPreparedViewportImage(
    SoBRLViewportImage &target, int layer, SbBool visible,
    BObolViewportGeometryPublication geometryPublication) :
    BObolPreparedViewportImage(target)
{
    this->next().layer = layer;
    this->next().visible = visible;
    (void)this->prepare(geometryPublication);
}

SoBRLViewportImage &
BObolPreparedViewportImage::next()
{
    if (this->impl->ready || this->impl->committed)
	throw std::logic_error("viewport presentation already prepared");
    return *static_cast<SoBRLViewportImage *>(
	this->impl->candidateOwner.get());
}

bool
BObolPreparedViewportImage::prepare(
    BObolViewportGeometryPublication geometryPublication)
{
    auto *candidate =
	static_cast<SoBRLViewportImage *>(this->impl->candidateOwner.get());
    if (geometryPublication == BObolViewportGeometryPublication::Preserve)
	return this->preparePayload(nullptr, false);

    struct bobol_image_payload payload;
    if (candidate->visible.getValue()) {
	const int loaded = geometryPublication ==
	    BObolViewportGeometryPublication::RebuildFromRetainedTexture ?
	    bobol_image_payload_load_retained(&this->impl->target, &payload) :
	    bobol_image_payload_load_current(
		this->impl->target.getImageSource(),
		&payload);
	if (loaded != 0)
	    return false;
    }
    return this->preparePayload(candidate->visible.getValue() ? &payload :
	nullptr, true);
}

bool
BObolPreparedViewportImage::prepare(
    const struct bobol_image_payload &payload)
{
    return this->preparePayload(&payload, true);
}

bool
BObolPreparedViewportImage::preparePayload(
    const struct bobol_image_payload *payload, bool publishGeometry)
{
    if (this->impl->ready || this->impl->committed)
	throw std::logic_error("viewport presentation already prepared");
    auto *candidate =
	static_cast<SoBRLViewportImage *>(this->impl->candidateOwner.get());
    if (publishGeometry) {
	candidate->imageSource.setValue(
	    this->impl->target.imageSource.getValue());
	if (candidate->visible.getValue()) {
	    if (!payload)
		return false;
	    SoTexture2 *texture = NULL;
	    SoFaceSet *face = NULL;
	    auto childOwner = bobol_viewport_image_make_geometry(*candidate,
		*payload, &texture, &face);
	    if (!childOwner)
		return false;
	    candidate->addChild(childOwner.get());
	    candidate->textureNode = texture;
	    candidate->imageFaceSet = face;
	    candidate->realizedDataRevision = payload->dataRevision;
	    candidate->realizedDirtyRevision = payload->dirtyRevision;
	}

	std::vector<SoNode *> children;
	children.reserve(static_cast<size_t>(candidate->getNumChildren()));
	for (int i = 0; i < candidate->getNumChildren(); ++i)
	    children.push_back(candidate->getChild(i));
	this->impl->children =
	    this->impl->target.getChildren()->prepareReplacement(children);
	this->impl->texture = candidate->getTextureNode();
	this->impl->face = candidate->getImageFaceSet();
    }
    this->impl->fields =
	std::make_unique<PreparedScalarFields>(this->impl->target, *candidate);
    this->impl->ready = true;
    return true;
}

BObolPreparedViewportImage::~BObolPreparedViewportImage() = default;

bool
BObolPreparedViewportImage::valid() const
{
    return this->impl->ready;
}

void
BObolPreparedViewportImage::commit()
{
    if (!this->impl->ready || this->impl->committed)
	throw std::logic_error("invalid viewport presentation commit");
    this->impl->fields->commit();
    if (this->impl->children) {
	this->impl->children->commit();
	this->impl->target.textureNode = this->impl->texture;
	this->impl->target.imageFaceSet = this->impl->face;
    }
    this->impl->committed = true;
}

void
BObolPreparedViewportImage::restore()
{
    if (this->impl->fields)
	this->impl->fields->restore();
}

void
BObolPreparedViewportImage::notify(std::exception_ptr &failure)
{
    if (!this->impl->committed)
	return;
    this->restore();
    if (this->impl->children) {
	try { this->impl->children->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
    this->impl->fields->notify(failure);
}

static constexpr size_t framebuffer_layer_root_count = 3;

struct BObolPreparedFramebufferRoots::Impl {
    std::vector<std::unique_ptr<SoChildList::Replacement>> replacements;
};

BObolPreparedFramebufferRoots::BObolPreparedFramebufferRoots(
    BObolViewController &controller, SoNode *viewport, SoGroup *destination,
    BObolFramebufferRootInsertion insertion) :
    impl(new Impl)
{
    const std::array<SoGroup *, framebuffer_layer_root_count> roots{{
	controller.getFramebufferUnderlayRoot(),
	controller.getFramebufferInterlayRoot(),
	controller.getFramebufferOverlayRoot()}};
    this->impl->replacements.reserve(roots.size());
    for (SoGroup *root : roots) {
	if (!root)
	    continue;
	std::vector<SoNode *> next;
	next.reserve(static_cast<size_t>(root->getNumChildren()) +
	    (root == destination ? 1u : 0u));
	if (root == destination &&
	    insertion == BObolFramebufferRootInsertion::First)
	    next.push_back(viewport);
	for (int i = 0; i < root->getNumChildren(); ++i) {
	    SoNode *child = root->getChild(i);
	    if (child != viewport)
		next.push_back(child);
	}
	if (root == destination &&
	    insertion == BObolFramebufferRootInsertion::Last)
	    next.push_back(viewport);

	bool unchanged =
	    next.size() == static_cast<size_t>(root->getNumChildren());
	for (size_t i = 0; unchanged && i < next.size(); ++i)
	    unchanged = root->getChild(static_cast<int>(i)) == next[i];
	if (!unchanged)
	    this->impl->replacements.push_back(
		root->getChildren()->prepareReplacement(next));
    }
}

BObolPreparedFramebufferRoots::~BObolPreparedFramebufferRoots() = default;

void
BObolPreparedFramebufferRoots::commit()
{
    for (auto &replacement : this->impl->replacements)
	replacement->commit();
}

void
BObolPreparedFramebufferRoots::notify(std::exception_ptr &failure)
{
    for (auto &replacement : this->impl->replacements) {
	try { replacement->notify(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
}

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

static int
image_payload_load(SoBRLImageSource *source, struct bobol_image_payload *payload,
    bool refreshSource)
{
    if (!source || !payload)
	return -1;

    if (refreshSource && source->refreshFromStream() != 0)
	return -1;

    imgstream_t *stream = source->getStream();
    if (!stream)
	return -1;

    struct imgstream_info info;
    if (imgstream_get_info(stream, &info) != 0)
	return -1;
    return bobol_image_payload_load_info(stream, info,
	image_display_clamp_u32(source->dirtyRevision.getValue()), payload);
}

int
bobol_image_payload_load_info(imgstream_t *stream,
    const struct imgstream_info &info, uint32_t dirtyRevision,
    struct bobol_image_payload *payload)
{
    if (!stream || !payload)
	return -1;
    if (info.width == 0 || info.height == 0 || info.channels == 0)
	return -1;
    if (info.width > (size_t)INT_MAX || info.height > (size_t)INT_MAX ||
	info.channels > (size_t)INT_MAX)
	return -1;
    if (info.width > SIZE_MAX / info.channels)
	return -1;

    const size_t stride = info.width * info.channels;
    if (info.height > SIZE_MAX / stride)
	return -1;

    payload->pixels.assign(stride * info.height, 0);
    if (imgstream_read_rect(stream, 0, 0, info.width, info.height,
	    payload->pixels.data(), stride) != 0) {
	payload->pixels.clear();
	return -1;
    }

    struct imgstream_info realized;
    if (imgstream_get_info(stream, &realized) != 0 ||
	realized.width != info.width || realized.height != info.height ||
	realized.format != info.format || realized.channels != info.channels ||
	realized.generation != info.generation || realized.dirty != info.dirty ||
	realized.producer_active != info.producer_active ||
	(info.dirty && (realized.dirty_rect.x != info.dirty_rect.x ||
	    realized.dirty_rect.y != info.dirty_rect.y ||
	    realized.dirty_rect.width != info.dirty_rect.width ||
	    realized.dirty_rect.height != info.dirty_rect.height))) {
	payload->pixels.clear();
	return -1;
    }

    payload->width = (int)info.width;
    payload->height = (int)info.height;
    payload->channels = (int)info.channels;
    payload->dataRevision = image_display_clamp_u32(info.generation);
    payload->dirtyRevision = dirtyRevision;
    return 0;
}

int
bobol_image_payload_load(SoBRLImageSource *source,
    struct bobol_image_payload *payload)
{
    return image_payload_load(source, payload, true);
}

int
bobol_image_payload_load_current(SoBRLImageSource *source,
    struct bobol_image_payload *payload)
{
    return image_payload_load(source, payload, false);
}

int
bobol_image_payload_load_retained(const SoBRLViewportImage *viewport,
    struct bobol_image_payload *payload)
{
    if (!viewport || !payload)
	return -1;
    SoTexture2 *texture = viewport->getTextureNode();
    int width = 0;
    int height = 0;
    int channels = 0;
    const unsigned char *pixels = texture ?
	texture->getImageData(width, height, channels) : nullptr;
    if (!pixels || width <= 0 || height <= 0 || channels <= 0)
	return -1;

    const size_t pixelWidth = static_cast<size_t>(width);
    const size_t pixelHeight = static_cast<size_t>(height);
    const size_t componentCount = static_cast<size_t>(channels);
    if (pixelWidth > SIZE_MAX / pixelHeight ||
	pixelWidth * pixelHeight > SIZE_MAX / componentCount)
	return -1;

    payload->width = width;
    payload->height = height;
    payload->channels = channels;
    payload->dataRevision = viewport->realizedDataRevision.getValue();
    payload->dirtyRevision = viewport->realizedDirtyRevision.getValue();
    const size_t pixelCount = pixelWidth * pixelHeight * componentCount;
    payload->pixels.assign(pixels, pixels + pixelCount);
    return 0;
}

void
bobol_image_fit_size(float sourceWidth, float sourceHeight,
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
bobol_image_texture_rect(float sourceWidth, float sourceHeight,
			   const SbVec2f &sourceCenter, float sourceZoom,
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

    if (u0)
	*u0 = outU0;
    if (v0)
	*v0 = outV0;
    if (u1)
	*u1 = outU0 + spanU;
    if (v1)
	*v1 = outV0 + spanV;
}

SbModernUtils::SoNodeRef
bobol_image_make_textured_quad(const struct bobol_image_payload *payload,
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
	return SbModernUtils::SoNodeRef(nullptr);

    SbModernUtils::SoNodeRef rootOwner(new SoSeparator);
    auto *root = static_cast<SoSeparator *>(rootOwner.get());

    SbModernUtils::SoNodeRef pickOwner(new SoPickStyle);
    auto *pick = static_cast<SoPickStyle *>(pickOwner.get());
    pick->style = selectable ? SoPickStyle::SHAPE : SoPickStyle::UNPICKABLE;
    root->addChild(pick);

    SbModernUtils::SoNodeRef depthOwner(new SoDepthBuffer);
    auto *depth = static_cast<SoDepthBuffer *>(depthOwner.get());
    depth->test = depthTest;
    depth->write = depthWrite;
    depth->function = SoDepthBuffer::LEQUAL;
    root->addChild(depth);

    SbModernUtils::SoNodeRef hintsOwner(new SoShapeHints);
    auto *hints = static_cast<SoShapeHints *>(hintsOwner.get());
    hints->vertexOrdering = SoShapeHints::COUNTERCLOCKWISE;
    hints->shapeType = doubleSided ? SoShapeHints::UNKNOWN_SHAPE_TYPE : SoShapeHints::SOLID;
    hints->faceType = SoShapeHints::CONVEX;
    root->addChild(hints);

    SbModernUtils::SoNodeRef textureOwner(new SoTexture2);
    auto *texture = static_cast<SoTexture2 *>(textureOwner.get());
    texture->model = SoTexture2::REPLACE;
    texture->wrapS = SoTexture2::CLAMP;
    texture->wrapT = SoTexture2::CLAMP;
    /* A new texture already has an empty filename. setImageData() writes that
     * value again, which can defer its filename sensor during observer reentry
     * until after the pixels are installed and then clear them as a reset. */
    texture->image.setValue(SbVec2s(payload->width, payload->height),
	payload->channels, payload->pixels.data());
    texture->image.setDefault(FALSE);
    root->addChild(texture);

    SbModernUtils::SoNodeRef materialOwner(new SoMaterial);
    auto *material = static_cast<SoMaterial *>(materialOwner.get());
    material->diffuseColor.setValue(1.0f, 1.0f, 1.0f);
    material->transparency.set1Value(0, 1.0f - image_display_clamp_float(opacity, 0.0f, 1.0f));
    root->addChild(material);

    SbModernUtils::SoNodeRef texCoordOwner(new SoTextureCoordinate2);
    auto *texCoord = static_cast<SoTextureCoordinate2 *>(texCoordOwner.get());
    texCoord->point.set1Value(0, SbVec2f(u0, v0));
    texCoord->point.set1Value(1, SbVec2f(u1, v0));
    texCoord->point.set1Value(2, SbVec2f(u1, v1));
    texCoord->point.set1Value(3, SbVec2f(u0, v1));
    root->addChild(texCoord);

    SbModernUtils::SoNodeRef coordsOwner(new SoCoordinate3);
    auto *coords = static_cast<SoCoordinate3 *>(coordsOwner.get());
    coords->point.set1Value(0, SbVec3f(x0, y0, z0));
    coords->point.set1Value(1, SbVec3f(x0 + width, y0, z0));
    coords->point.set1Value(2, SbVec3f(x0 + width, y0 + height, z0));
    coords->point.set1Value(3, SbVec3f(x0, y0 + height, z0));
    root->addChild(coords);

    SbModernUtils::SoNodeRef faceOwner(new SoFaceSet);
    auto *face = static_cast<SoFaceSet *>(faceOwner.get());
    face->numVertices.setValue(4);
    root->addChild(face);

    if (textureOut)
	*textureOut = texture;
    if (faceOut)
	*faceOut = face;
    return rootOwner;
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

SbModernUtils::SoNodeRef
bobol_viewport_image_make_geometry(const SoBRLViewportImage &viewport,
    const struct bobol_image_payload &payload, SoTexture2 **textureOut,
    SoFaceSet **faceOut)
{
    float displayWidth = 0.0f;
    float displayHeight = 0.0f;
    const SbVec2f requested = viewport.size.getValue();
    bobol_image_fit_size(static_cast<float>(payload.width),
	static_cast<float>(payload.height), requested[0], requested[1],
	viewport.fit.getValue(), viewport.preserveAspect.getValue() == TRUE,
	&displayWidth, &displayHeight);

    const SbVec2f position = viewport.position.getValue();
    float x0 = 0.0f;
    float y0 = 0.0f;
    viewport_anchor_origin(viewport.anchor.getValue(), position[0], position[1],
	displayWidth, displayHeight, &x0, &y0);

    float u0 = 0.0f;
    float v0 = 0.0f;
    float u1 = 1.0f;
    float v1 = 1.0f;
    bobol_image_texture_rect(static_cast<float>(payload.width),
	static_cast<float>(payload.height), viewport.sourceCenter.getValue(),
	viewport.sourceZoom.getValue(),
	&u0, &v0, &u1, &v1);

    SoTexture2 *texture = NULL;
    SoFaceSet *face = NULL;
    auto quadOwner = bobol_image_make_textured_quad(&payload, x0, y0,
	static_cast<float>(viewport.zOrder.getValue()), displayWidth,
	displayHeight, u0, v0, u1, v1, viewport.opacity.getValue(), FALSE,
	FALSE, FALSE, TRUE, &texture, &face);
    if (!quadOwner)
	return SbModernUtils::SoNodeRef(nullptr);

    SbModernUtils::SoNodeRef hudOwner(new SoHUDKit);
    static_cast<SoHUDKit *>(hudOwner.get())->addWidget(quadOwner.get());
    if (textureOut)
	*textureOut = texture;
    if (faceOut)
	*faceOut = face;
    return hudOwner;
}

void
bobol_image_publish_geometry(SoSeparator &target, SoNode *child,
	SoTexture2 *texture, SoFaceSet *face,
	SoTexture2 *&textureSlot, SoFaceSet *&faceSlot,
	SoSFUInt32 &dataRevision, SoSFUInt32 &dirtyRevision,
	uint32_t nextDataRevision, uint32_t nextDirtyRevision,
	bool publishRevisions)
{
    std::vector<SoNode *> children;
    if (child)
	children.push_back(child);
    auto replacement = target.getChildren()->prepareReplacement(children);
    if (!publishRevisions) {
	replacement->commit();
	textureSlot = NULL;
	faceSlot = NULL;
	replacement->notify();
	return;
    }

    PreparedFieldNotifications<2> notifications(target, {{
	{&dataRevision, true}, {&dirtyRevision, true}
    }});
    replacement->commit();
    textureSlot = texture;
    faceSlot = face;
    dataRevision = nextDataRevision;
    dirtyRevision = nextDirtyRevision;

    notifications.restore();
    std::exception_ptr failure;
    try { replacement->notify(); }
    catch (...) { failure = std::current_exception(); }
    notifications.notify(failure);
    if (failure)
	std::rethrow_exception(failure);
}
