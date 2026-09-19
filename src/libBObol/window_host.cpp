/*                W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BImageSource.h"
#include "BObol/BFramebuffer.h"
#include "BObol/BViewController.h"
#include "BObol/BViewportImage.h"
#include "BObol/BWindowHost.h"
#include "image_display_util.h"
#include "image_source_private.h"
#include "scalar_publication_private.h"
#include "view_controller_private.h"
#include "window_host_private.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <string.h>
#include <vector>

struct BObolFramebufferAttachment {
    imgstream_fb_t *fb;
    SoBRLImageSource *source;
    SoBRLViewportImage *viewport;
    BObolFramebufferComposition composition;
};

static constexpr unsigned int minimum_window_extent = 1;
static constexpr const char *default_window_title = "BRL-CAD Obol";

static bool
framebuffer_viewport_realizes_source(const SoBRLViewportImage &viewport,
    const SoBRLImageSource &source)
{
    if (!viewport.visible.getValue())
	return viewport.getNumChildren() == 0 && !viewport.getTextureNode() &&
	    !viewport.getImageFaceSet();
    return viewport.getNumChildren() == 1 && viewport.getTextureNode() &&
	viewport.getImageFaceSet() &&
	viewport.realizedDataRevision.getValue() ==
	    source.dataRevision.getValue() &&
	viewport.realizedDirtyRevision.getValue() ==
	    source.dirtyRevision.getValue();
}

static BObolWindowDesc default_desc(void);

struct BObolWindowHostPrivate {
    BObolWindowHostPrivate(void) :
	controller(NULL),
	ownsController(TRUE),
	open(FALSE),
	desc(std::make_unique<BObolWindowDesc>(default_desc())),
	pollRate(0)
    {
	auto initialController = std::make_unique<BObolViewController>();
	controller = initialController.release();
    }

    BObolViewController *controller;
    SbBool ownsController;
    SbBool open;
    /* Swapping the prepared descriptor keeps failed opens from exposing a
     * partially copied policy through getDesc(). */
    std::unique_ptr<BObolWindowDesc> desc;
    long pollRate;
    BObolInputContext input;
    std::vector<BObolFramebufferAttachment> framebuffers;
    std::vector<BObolFramebufferStream *> framebufferStreams;
    std::vector<imgstream_fb_t *> displayFramebuffers;
};

static BObolWindowDesc
default_desc(void)
{
    BObolWindowDesc desc;
    desc.mode = BOBOL_WINDOW_HEADLESS;
    desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = minimum_window_extent;
    desc.height = minimum_window_extent;
    desc.title = default_window_title;
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;
    return desc;
}

static unsigned int
normalized_window_extent(unsigned int extent)
{
    return extent ? extent : minimum_window_extent;
}

static void
sanitize_desc(BObolWindowDesc *desc)
{
    if (!desc)
	return;
    desc->width = normalized_window_extent(desc->width);
    desc->height = normalized_window_extent(desc->height);
}

static bool
desc_matches_request(const BObolWindowDesc &current,
	const BObolWindowDesc *requested)
{
    const BObolWindowMode mode = requested ? requested->mode :
	BOBOL_WINDOW_HEADLESS;
    const BObolWindowBackend backend = requested ? requested->backend :
	BOBOL_WINDOW_BACKEND_OFFSCREEN;
    const unsigned int width = requested ?
	normalized_window_extent(requested->width) : minimum_window_extent;
    const unsigned int height = requested ?
	normalized_window_extent(requested->height) : minimum_window_extent;
    const char *title = requested ? requested->title.getString() :
	default_window_title;
    const char *display = requested ? requested->display.getString() : "";
    const char *nativeIdHint = requested ?
	requested->nativeIdHint.getString() : "";
    const SbBool visible = requested ? requested->visible : FALSE;
    return current.mode == mode && current.backend == backend &&
	current.width == width && current.height == height &&
	current.title == title && current.display == display &&
	current.nativeIdHint == nativeIdHint && current.visible == visible;
}

static BObolWindowBackend
backend_from_fb_display(enum imgstream_fb_display_kind display)
{
    switch (display) {
	case IMGSTREAM_FB_DISPLAY_QTGL:
	    return BOBOL_WINDOW_BACKEND_QT;
	case IMGSTREAM_FB_DISPLAY_SWRAST:
	    return BOBOL_WINDOW_BACKEND_OFFSCREEN;
	case IMGSTREAM_FB_DISPLAY_X:
	    return BOBOL_WINDOW_BACKEND_TK;
	case IMGSTREAM_FB_DISPLAY_OGL:
	case IMGSTREAM_FB_DISPLAY_WGL:
	    return BOBOL_WINDOW_BACKEND_OPENGL;
	case IMGSTREAM_FB_DISPLAY_NONE:
	default:
	    return BOBOL_WINDOW_BACKEND_AUTO;
    }
}

SoGroup *
bobol_window_host_root_group(BObolViewController *controller)
{
    if (!controller)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<SoGroup *>(root);
}

static void
request_framebuffer_presentation(BObolViewController *controller,
	const char *reason)
{
    /* Framebuffer pixels and placement are presentation state.  They must be
     * displayed promptly, but they neither invalidate CAD demand nor provide
     * a renderer-capacity sample for the retained geometry allocator. */
    if (controller)
	controller->requestPresentationRender(reason);
}

static float
positive_zoom(int xzoom, int yzoom)
{
    int x = xzoom > 0 ? xzoom : 1;
    int y = yzoom > 0 ? yzoom : 1;
    return (float)((x < y) ? x : y);
}

static bool
same_float(float a, float b)
{
    return std::fabs(a - b) <= 1.0e-6f;
}

template <typename PrepareRequest, typename CommitRequest,
    typename NotifyRequest>
static int
publish_retained_framebuffer_viewport(
	BObolPreparedViewportImage &viewportPublication,
	PrepareRequest prepareRequest, CommitRequest commitRequest,
	NotifyRequest notifyRequest)
{
    /* Presentation inputs retain the displayed pixels and revisions. Flush
     * owns delivery of a pending stream generation. */
    if (!viewportPublication.prepare(
	    BObolViewportGeometryPublication::RebuildFromRetainedTexture))
	return -1;
    prepareRequest();

    viewportPublication.commit();
    commitRequest();

    std::exception_ptr failure;
    viewportPublication.notify(failure);
    try { notifyRequest(); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
    return 0;
}

static BObolFramebufferAttachment *
find_attachment(std::vector<BObolFramebufferAttachment> &attachments,
		imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

static const BObolFramebufferAttachment *
find_attachment_const(const std::vector<BObolFramebufferAttachment> &attachments,
		      imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

BObolWindowHost::BObolWindowHost(void) :
    p(new BObolWindowHostPrivate)
{
}

BObolWindowHost::~BObolWindowHost(void)
{
    try { this->detachDisplayFramebuffers(); }
    catch (...) {}
    try { this->detachFramebufferStreams(); }
    catch (...) {}
    this->destroyFramebuffersNoexcept();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = NULL;
    delete this->p;
}

void
BObolWindowHost::registerFramebufferStream(BObolFramebufferStream *stream)
{
    if (!stream || std::find(this->p->framebufferStreams.begin(),
	this->p->framebufferStreams.end(), stream) !=
	this->p->framebufferStreams.end())
	return;
    this->p->framebufferStreams.push_back(stream);
}

void
BObolWindowHost::unregisterFramebufferStream(BObolFramebufferStream *stream)
{
    this->p->framebufferStreams.erase(std::remove_if(
	this->p->framebufferStreams.begin(), this->p->framebufferStreams.end(),
	[stream](BObolFramebufferStream *candidate) {
	    return candidate == stream;
	}),
	this->p->framebufferStreams.end());
}

void
BObolWindowHost::registerDisplayFramebuffer(imgstream_fb_t *fb)
{
    if (!fb || std::find(this->p->displayFramebuffers.begin(),
	this->p->displayFramebuffers.end(), fb) != this->p->displayFramebuffers.end())
	return;
    this->p->displayFramebuffers.push_back(fb);
}

void
BObolWindowHost::unregisterDisplayFramebuffer(imgstream_fb_t *fb)
{
    this->p->displayFramebuffers.erase(std::remove_if(
	this->p->displayFramebuffers.begin(), this->p->displayFramebuffers.end(),
	[fb](imgstream_fb_t *candidate) {
	    return candidate == fb;
	}),
	this->p->displayFramebuffers.end());
}

void
BObolWindowHost::detachFramebufferStreams(void)
{
    for (size_t i = 0; i < this->p->framebufferStreams.size(); i++)
	this->p->framebufferStreams[i]->detachHost(this);
    this->p->framebufferStreams.clear();
}

void
BObolWindowHost::detachDisplayFramebuffers(void)
{
    while (!this->p->displayFramebuffers.empty()) {
	imgstream_fb_t *fb = this->p->displayFramebuffers.back();
	this->p->displayFramebuffers.pop_back();
	(void)imgstream_fb_detach_display_host(fb);
    }
}

int
BObolWindowHost::open(const BObolWindowDesc *desc)
{
    PreparedOpenPublication publication;
    if (this->prepareOpenPublication(desc, publication) != 0)
	return -1;
    this->commitOpenPublication(publication);
    this->notifyOpenPublication(publication);
    return 0;
}

int
BObolWindowHost::prepareOpenPublication(const BObolWindowDesc *desc,
	PreparedOpenPublication &publication)
{
    if (this->p->open && desc_matches_request(*this->p->desc, desc))
	return 0;
    auto nextDesc = std::make_unique<BObolWindowDesc>(
	desc ? *desc : default_desc());
    sanitize_desc(nextDesc.get());

    if (!this->p->controller ||
	!bobol_window_host_root_group(this->p->controller))
	return -1;

    this->p->controller->prepareViewportSizePublication(nextDesc->width,
	nextDesc->height, "viewport-size", publication.viewport);
    publication.desc = std::move(nextDesc);
    publication.changed = true;
    return 0;
}

void
BObolWindowHost::commitOpenPublication(
    PreparedOpenPublication &publication) noexcept
{
    if (!publication.changed)
	return;
    this->p->desc.swap(publication.desc);
    this->p->open = TRUE;
    this->p->controller->commitViewportRegionPublication(
	publication.viewport);
}

void
BObolWindowHost::notifyOpenPublication(
    const PreparedOpenPublication &publication)
{
    if (!publication.changed)
	return;
    this->p->controller->notifyViewportRegionPublication(
	publication.viewport);
}

SbBool
BObolWindowHost::openPublicationChangesViewport(
    const PreparedOpenPublication &publication) const
{
    return publication.viewport.changed ? TRUE : FALSE;
}

void
BObolWindowHost::close(void)
{
    this->prepareClosePublication();
    this->commitClosePublication();
}

void
BObolWindowHost::prepareClosePublication(void)
{
    while (!this->p->framebuffers.empty())
	this->closeFramebuffer(this->p->framebuffers.back().fb);
}

void
BObolWindowHost::commitClosePublication(void) noexcept
{
    this->p->open = FALSE;
}

void
BObolWindowHost::destroyFramebuffersNoexcept(void) noexcept
{
    while (!this->p->framebuffers.empty()) {
	const size_t precedingCount = this->p->framebuffers.size();
	try {
	    this->closeFramebuffer(this->p->framebuffers.back().fb);
	} catch (...) {
	}
	if (this->p->framebuffers.size() == precedingCount)
	    break;
    }

    /* Preparation failure leaves a complete attachment in its roots. Release
     * the host's manual references; an owning controller retires the remaining
     * root references below, while a borrowed controller keeps valid nodes. */
    for (const BObolFramebufferAttachment &attachment : this->p->framebuffers) {
	try { attachment.viewport->unref(); }
	catch (...) {}
	try { attachment.source->unref(); }
	catch (...) {}
    }
    this->p->framebuffers.clear();
    this->p->open = FALSE;
}

SbBool
BObolWindowHost::isOpen(void) const
{
    return this->p->open;
}

const BObolWindowDesc &
BObolWindowHost::getDesc(void) const
{
    return *this->p->desc;
}

SbBool
BObolWindowHost::attachController(BObolViewController *controller,
				    SbBool takeOwnership)
{
    if (controller && !bobol_window_host_root_group(controller))
	return FALSE;
    if (controller == this->p->controller) {
	this->p->ownsController = controller && takeOwnership ? TRUE : FALSE;
	return TRUE;
    }

    this->close();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = controller;
    this->p->ownsController = controller && takeOwnership ? TRUE : FALSE;
    return TRUE;
}

void
bobol_window_host_detach_controller_noexcept(BObolWindowHost *host) noexcept
{
    if (!host || !host->p)
	return;
    try { host->close(); }
    catch (...) {}
    host->destroyFramebuffersNoexcept();
    /* Display endpoints always borrow their controller into a host. Terminal
     * endpoint cleanup must sever that pointer even when a derived close path
     * failed, without invoking another live publication. */
    host->p->controller = NULL;
    host->p->ownsController = FALSE;
}

BObolViewController *
BObolWindowHost::getController(void) const
{
    return this->p->controller;
}

static int
window_host_input_action(void *data, BObolInputAction action,
	const BObolInputEvent *event)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->applyInputAction(action, event) : -1;
}

int
BObolWindowHost::handleInputEvent(const BObolInputEvent *event,
				    const BObolInputProfile *profile)
{
    if (!event)
	return -1;
    if (!profile || !profile->bindings || profile->bindingCount == 0)
	return 0;

    this->p->input.setProfile(profile);
    this->p->input.setActionHandler(window_host_input_action, this);
    return this->p->input.dispatch(event);
}

int
BObolWindowHost::applyInputAction(BObolInputAction action,
				    const BObolInputEvent *UNUSED(event))
{
    if (action == BOBOL_ACTION_NONE)
	return 0;
    if (!this->p->controller)
	return -1;

    this->p->controller->requestLodCapacityRender("input-action");
    return 1;
}

int
BObolWindowHost::poll(const BObolInputProfile *UNUSED(profile))
{
    return 0;
}

long
BObolWindowHost::pollRate(void) const
{
    return this->p->pollRate;
}

int
BObolWindowHost::openFramebuffer(imgstream_fb_t *fb,
				   const imgstream_fb_spec_info_t *info)
{
    if (!fb)
	return -1;
    if (find_attachment(this->p->framebuffers, fb))
	return 0;

    BObolWindowDesc desc = this->p->open ? *this->p->desc : default_desc();
    desc.width = (unsigned int)std::max<size_t>(imgstream_fb_width(fb), 1);
    desc.height = (unsigned int)std::max<size_t>(imgstream_fb_height(fb), 1);
    desc.mode = (info && info->display == IMGSTREAM_FB_DISPLAY_SWRAST) ?
		BOBOL_WINDOW_HEADLESS : BOBOL_WINDOW_TOPLEVEL;
    desc.backend = info ? backend_from_fb_display(info->display) :
		   BOBOL_WINDOW_BACKEND_AUTO;
    desc.visible = desc.mode == BOBOL_WINDOW_TOPLEVEL;
    if (info && info->host && info->host_len)
	desc.display = SbString(info->host, 0, (int)info->host_len - 1);
    desc.title = imgstream_fb_name(fb) ? imgstream_fb_name(fb) : "BRL-CAD framebuffer";

    imgstream_t *stream = imgstream_fb_stream(fb);
    if (!stream)
	return -1;

    if (!bobol_window_host_root_group(this->p->controller) ||
	!this->p->controller->getFramebufferOverlayRoot())
	return -1;

    SbModernUtils::SoNodeRef sourceOwner(new SoBRLImageSource);
    auto *source = static_cast<SoBRLImageSource *>(sourceOwner.get());
    source->imageId = imgstream_fb_name(fb) ? imgstream_fb_name(fb) : "framebuffer";
    source->sourceUri = source->imageId.getValue();
    if (source->setStream(stream) != 0)
	return -1;

    SbModernUtils::SoNodeRef viewportOwner(new SoBRLViewportImage);
    auto *viewport = static_cast<SoBRLViewportImage *>(viewportOwner.get());
    viewport->overlayId = source->imageId.getValue();
    viewport->imageSource.setValue(source);
    viewport->layer = SoBRLViewportImage::OVERLAY;
    viewport->anchor = SoBRLViewportImage::LOWER_LEFT;
    viewport->fit = SoBRLViewportImage::STRETCH;
    viewport->preserveAspect = FALSE;
    viewport->position.setValue(0.0f, 0.0f);
    viewport->size.setValue((float)desc.width, (float)desc.height);
    viewport->sourceCenter.setValue((float)desc.width * 0.5f,
				    (float)desc.height * 0.5f);
    viewport->sourceZoom = 1.0f;
    viewport->cursorVisible = FALSE;
    if (viewport->rebuildGeometry() != 0)
	return -1;

    BObolFramebufferAttachment attachment;
    attachment.fb = fb;
    attachment.source = source;
    attachment.viewport = viewport;
    attachment.composition = BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY;
    if (this->p->framebuffers.size() == this->p->framebuffers.max_size())
	return -1;
    this->p->framebuffers.reserve(this->p->framebuffers.size() + 1);
    BObolPreparedFramebufferRoots roots(*this->p->controller, viewport,
	this->p->controller->getFramebufferOverlayRoot(),
	BObolFramebufferRootInsertion::Last);

    if (this->open(&desc) != 0)
	return -1;

    this->p->framebuffers.push_back(attachment);
    roots.commit();
    sourceOwner.release();
    viewportOwner.release();

    std::exception_ptr failure;
    roots.notify(failure);
    try { request_framebuffer_presentation(this->p->controller, "fb-open"); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
    return 0;
}

void
BObolWindowHost::closeFramebuffer(imgstream_fb_t *fb)
{
    for (size_t i = 0; i < this->p->framebuffers.size(); i++) {
	if (this->p->framebuffers[i].fb != fb)
	    continue;

	const BObolFramebufferAttachment closing = this->p->framebuffers[i];
	BObolPreparedFramebufferRoots roots(*this->p->controller,
	    closing.viewport, nullptr, BObolFramebufferRootInsertion::Last);
	this->p->framebuffers.erase(this->p->framebuffers.begin() + (ptrdiff_t)i);
	roots.commit();

	std::exception_ptr failure;
	roots.notify(failure);
	try { closing.viewport->unref(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
	try { closing.source->unref(); }
	catch (...) { if (!failure) failure = std::current_exception(); }
	try { request_framebuffer_presentation(this->p->controller, "fb-close"); }
	catch (...) { if (!failure) failure = std::current_exception(); }
	if (failure)
	    std::rethrow_exception(failure);
	return;
    }
}

int
BObolWindowHost::setFramebufferComposition(imgstream_fb_t *fb,
	BObolFramebufferComposition composition)
{
    if (composition < BOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	composition > BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	return -1;

    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    BObolViewController *controller = this->p->controller;
    if (!attachment || !controller)
	return -1;
    /*
     * The framebuffer bridge may synchronize its policy on every canvas
     * presentation.  Reapplying an unchanged mode used to remove, rebuild,
     * and reinsert the viewport and leave another "fb-composition" render
     * request behind each OSMesa frame.  Besides wasting owner-thread work,
     * that makes a fully settled progressive view observably never idle.
     */
    if (attachment->composition == composition)
	return 0;

    SoBRLViewportImage *viewport = attachment->viewport;
    SoGroup *destination = NULL;
    int viewportLayer = viewport->layer.getValue();
    SbBool visible = TRUE;
    switch (composition) {
	case BOBOL_FRAMEBUFFER_COMPOSITION_OFF:
	    visible = FALSE;
	    break;
	case BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY:
	    destination = controller->getFramebufferUnderlayRoot();
	    viewportLayer = SoBRLViewportImage::UNDERLAY;
	    break;
	case BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY:
	    destination = controller->getFramebufferInterlayRoot();
	    viewportLayer = SoBRLViewportImage::INTERLAY;
	    break;
	case BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY:
	    destination = controller->getFramebufferOverlayRoot();
	    viewportLayer = SoBRLViewportImage::OVERLAY;
	    break;
	default:
	    return -1;
    }
    if (composition != BOBOL_FRAMEBUFFER_COMPOSITION_OFF && !destination)
	return -1;

    BObolPreparedViewportImage viewportPublication(*viewport,
	viewportLayer, visible);
    if (!viewportPublication.valid())
	return -1;
    BObolPreparedFramebufferRoots roots(*controller, viewport, destination,
	BObolFramebufferRootInsertion::Last);

    viewportPublication.commit();
    roots.commit();
    attachment->composition = composition;

    std::exception_ptr failure;
    viewportPublication.notify(failure);
    roots.notify(failure);
    try { request_framebuffer_presentation(controller, "fb-composition"); }
    catch (...) { if (!failure) failure = std::current_exception(); }
    if (failure)
	std::rethrow_exception(failure);
    return 0;
}

int
BObolWindowHost::flushFramebuffer(imgstream_fb_t *fb)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    SoBRLImageSource *source = attachment->source;
    SoBRLViewportImage *viewport = attachment->viewport;
    BObolImageSourceRefresh sourceRefresh =
	BObolImageSourcePublication::refreshRequired(*source);
    if (sourceRefresh == BObolImageSourceRefresh::Failed)
	return -1;

    const bool retainedCurrent =
	framebuffer_viewport_realizes_source(*viewport, *source);
    if (sourceRefresh == BObolImageSourceRefresh::Current && retainedCurrent)
	return 0;

    std::unique_ptr<BObolImageSourcePublication> sourcePublication;
    if (sourceRefresh == BObolImageSourceRefresh::Required) {
	sourcePublication =
	    std::make_unique<BObolImageSourcePublication>(*source);
	sourceRefresh = sourcePublication->prepareRefresh();
	if (sourceRefresh == BObolImageSourceRefresh::Failed)
	    return -1;
    }

    const bool sourceChanged =
	sourceRefresh == BObolImageSourceRefresh::Required;
    const SoBRLImageSource &nextSource = sourceChanged ?
	sourcePublication->successor() : *source;
    const bool rebuildViewport =
	!framebuffer_viewport_realizes_source(*viewport, nextSource);
    if (!sourceChanged && !rebuildViewport)
	return 0;

    std::unique_ptr<BObolPreparedViewportImage> viewportPublication;
    if (rebuildViewport) {
	struct bobol_image_payload payload;
	if (viewport->visible.getValue()) {
	    const int loaded = sourceChanged ?
		sourcePublication->loadPreparedPayload(payload) :
		bobol_image_payload_load_current(source, &payload);
	    if (loaded != 0)
		return -1;
	}
	viewportPublication =
	    std::make_unique<BObolPreparedViewportImage>(*viewport);
	if (!viewportPublication->prepare(payload))
	    return -1;
    }

    BObolViewController *controller = this->p->controller;
    BObolPreparedRenderRequest renderRequest;
    if (controller)
	renderRequest = controller->prepareRenderRequest("fb-flush",
	    BObolViewController::RenderRequestIntent::PRESENTATION);

    if (sourceChanged)
	sourcePublication->commit();
    if (viewportPublication)
	viewportPublication->commit();
    if (controller)
	controller->commitRenderRequest(renderRequest);

    if (sourceChanged)
	sourcePublication->restore();
    if (viewportPublication)
	viewportPublication->restore();

    std::exception_ptr failure;
    if (sourceChanged)
	sourcePublication->notify(failure);
    if (viewportPublication)
	viewportPublication->notify(failure);
    try {
	if (controller)
	    controller->notifyRenderRequest(renderRequest);
    } catch (...) {
	if (!failure)
	    failure = std::current_exception();
    }
    if (failure)
	std::rethrow_exception(failure);
    return 0;
}

int
BObolWindowHost::resetFramebuffer(imgstream_fb_t *fb)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    const SbVec2f center(
	static_cast<float>(imgstream_fb_width(fb)) * 0.5f,
	static_cast<float>(imgstream_fb_height(fb)) * 0.5f);
    SoBRLViewportImage *viewport = attachment->viewport;
    if (viewport->sourceCenter.getValue() == center &&
	same_float(viewport->sourceZoom.getValue(), 1.0f) &&
	viewport->cursorVisible.getValue() == FALSE)
	return 0;

    BObolPreparedViewportImage viewportPublication(*viewport);
    viewportPublication.next().sourceCenter = center;
    viewportPublication.next().sourceZoom = 1.0f;
    viewportPublication.next().cursorVisible = FALSE;
    BObolViewController *controller = this->p->controller;
    BObolPreparedRenderRequest renderRequest;
    return publish_retained_framebuffer_viewport(viewportPublication,
	[&] {
	    if (controller)
		renderRequest = controller->prepareRenderRequest("fb-reset",
		    BObolViewController::RenderRequestIntent::PRESENTATION);
	},
	[&] {
	    if (controller)
		controller->commitRenderRequest(renderRequest);
	},
	[&] {
	    if (controller)
		controller->notifyRenderRequest(renderRequest);
	});
}

int
BObolWindowHost::setFramebufferViewport(imgstream_fb_t *fb,
	int left, int top, int right, int bottom)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    int width = right - left;
    int height = bottom - top;
    if (width <= 0)
	width = (int)imgstream_fb_width(fb);
    if (height <= 0)
	height = (int)imgstream_fb_height(fb);
    if (width <= 0)
	width = 1;
    if (height <= 0)
	height = 1;

    SoBRLViewportImage *viewport = attachment->viewport;
    const SbVec2f position(static_cast<float>(left),
	static_cast<float>(top));
    const SbVec2f size(static_cast<float>(width),
	static_cast<float>(height));
    const bool viewportChanged = viewport->position.getValue() != position ||
	viewport->size.getValue() != size;
    BObolViewController *controller = this->p->controller;
    bool controllerChanged = false;
    if (controller) {
	const SbVec2s controllerSize =
	    controller->getViewportRegion().getWindowSize();
	controllerChanged = controllerSize[0] != width ||
	    controllerSize[1] != height;
    }
    if (!viewportChanged && !controllerChanged)
	return 0;

    if (!viewportChanged) {
	BObolViewController::PreparedViewportPublication publication;
	controller->prepareViewportSizePublication(
	    static_cast<unsigned int>(width), static_cast<unsigned int>(height),
	    "fb-viewport", publication);
	controller->commitViewportRegionPublication(publication);
	controller->notifyViewportRegionPublication(publication);
	return 0;
    }

    BObolPreparedViewportImage viewportPublication(*viewport);
    viewportPublication.next().position = position;
    viewportPublication.next().size = size;
    BObolViewController::PreparedViewportPublication controllerPublication;
    BObolPreparedRenderRequest placementRequest;
    return publish_retained_framebuffer_viewport(viewportPublication,
	[&] {
	    if (!controller)
		return;
	    controller->prepareViewportSizePublication(
		static_cast<unsigned int>(width),
		static_cast<unsigned int>(height), "fb-viewport",
		controllerPublication);
	    if (!controllerPublication.changed)
		placementRequest = controller->prepareRenderRequest("fb-viewport",
		    BObolViewController::RenderRequestIntent::PRESENTATION);
	},
	[&] {
	    if (!controller)
		return;
	    if (controllerPublication.changed)
		controller->commitViewportRegionPublication(
		    controllerPublication);
	    else
		controller->commitRenderRequest(placementRequest);
	},
	[&] {
	    if (!controller)
		return;
	    if (controllerPublication.changed)
		controller->notifyViewportRegionPublication(
		    controllerPublication);
	    else
		controller->notifyRenderRequest(placementRequest);
	});
}

int
BObolWindowHost::setFramebufferView(imgstream_fb_t *fb,
				      const imgstream_fb_view_t *view)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment || !view)
	return -1;

    float zoom = positive_zoom(view->xzoom, view->yzoom);
    SbVec2f oldCenter = attachment->viewport->sourceCenter.getValue();
    if (same_float(oldCenter[0], (float)view->xcenter) &&
	    same_float(oldCenter[1], (float)view->ycenter) &&
	    same_float(attachment->viewport->sourceZoom.getValue(), zoom))
	return 0;

    BObolPreparedViewportImage viewportPublication(*attachment->viewport);
    viewportPublication.next().sourceCenter.setValue(
	static_cast<float>(view->xcenter), static_cast<float>(view->ycenter));
    viewportPublication.next().sourceZoom = zoom;
    BObolViewController *controller = this->p->controller;
    BObolPreparedRenderRequest renderRequest;
    return publish_retained_framebuffer_viewport(viewportPublication,
	[&] {
	    if (controller)
		renderRequest = controller->prepareRenderRequest("fb-view",
		    BObolViewController::RenderRequestIntent::PRESENTATION);
	},
	[&] {
	    if (controller)
		controller->commitRenderRequest(renderRequest);
	},
	[&] {
	    if (controller)
		controller->notifyRenderRequest(renderRequest);
	});
}

int
BObolWindowHost::publishFramebufferCursorState(SoBRLViewportImage *viewport,
	SbBool visible, float x, float y, int shape, const char *reason)
{
    if (!viewport)
	return -1;

    const SbVec2f position(x, y);
    if (viewport->cursorVisible.getValue() == visible &&
	viewport->cursorImagePosition.getValue() == position &&
	viewport->cursorShape.getValue() == shape)
	return 0;

    BObolViewController *controller = this->p->controller;
    BObolPreparedRenderRequest renderRequest = controller ?
	controller->prepareRenderRequest(reason,
	    BObolViewController::RenderRequestIntent::PRESENTATION) :
	BObolPreparedRenderRequest();
    PreparedFieldNotifications<3> fields(*viewport, {{
	{&viewport->cursorVisible,
	    viewport->cursorVisible.getValue() != visible},
	{&viewport->cursorImagePosition,
	    (viewport->cursorImagePosition.getValue() != position) != FALSE},
	{&viewport->cursorShape,
	    viewport->cursorShape.getValue() != shape}
    }});

    viewport->cursorVisible = visible;
    viewport->cursorImagePosition = position;
    viewport->cursorShape = shape;
    if (controller)
	controller->commitRenderRequest(renderRequest);

    fields.restore();
    std::exception_ptr failure;
    fields.notify(failure);
    if (controller) {
	try { controller->notifyRenderRequest(renderRequest); }
	catch (...) { if (!failure) failure = std::current_exception(); }
    }
    if (failure)
	std::rethrow_exception(failure);
    return 0;
}

int
BObolWindowHost::setFramebufferCursor(imgstream_fb_t *fb,
					const imgstream_fb_cursor_t *cursor)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment || !cursor)
	return -1;

    return this->publishFramebufferCursorState(attachment->viewport,
	cursor->mode ? TRUE : FALSE, static_cast<float>(cursor->x),
	static_cast<float>(cursor->y), cursor->mode ?
	    SoBRLViewportImage::CURSOR_DEFAULT : SoBRLViewportImage::CURSOR_NONE,
	"fb-cursor");
}

int
BObolWindowHost::setFramebufferScreenCursor(imgstream_fb_t *fb,
	int mode, int x, int y)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    return this->publishFramebufferCursorState(attachment->viewport,
	mode ? TRUE : FALSE, static_cast<float>(x), static_cast<float>(y),
	mode ? SoBRLViewportImage::CURSOR_DEFAULT :
	    SoBRLViewportImage::CURSOR_NONE, "fb-screen-cursor");
}

int
BObolWindowHost::setFramebufferCursorShape(imgstream_fb_t *fb,
	const unsigned char *bits, int xbits, int ybits,
	int UNUSED(xorig), int UNUSED(yorig))
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    SbBool visible = attachment->viewport->cursorVisible.getValue();
    int shape = visible ? SoBRLViewportImage::CURSOR_DEFAULT :
	SoBRLViewportImage::CURSOR_NONE;
    if (bits && xbits > 0 && ybits > 0) {
	visible = TRUE;
	shape = SoBRLViewportImage::CURSOR_CUSTOM;
    }
    const SbVec2f position =
	attachment->viewport->cursorImagePosition.getValue();
    return this->publishFramebufferCursorState(attachment->viewport, visible,
	position[0], position[1], shape, "fb-cursor-shape");
}

int
BObolWindowHost::getFramebufferCount(void) const
{
    return (int)this->p->framebuffers.size();
}

SoBRLImageSource *
BObolWindowHost::getFramebufferImageSource(imgstream_fb_t *fb) const
{
    const BObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->source : NULL;
}

SoBRLViewportImage *
BObolWindowHost::getFramebufferViewportImage(imgstream_fb_t *fb) const
{
    const BObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->viewport : NULL;
}

struct BObolFramebufferBridge {
    static int
    open(imgstream_fb_t *fb, const imgstream_fb_spec_info_t *info, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    if (!host || host->openFramebuffer(fb, info) != 0)
	return -1;
    host->registerDisplayFramebuffer(fb);
    return 0;
}

    static void
    close(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    if (host) {
	host->unregisterDisplayFramebuffer(fb);
	host->closeFramebuffer(fb);
    }
}

    static int
    flush(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->flushFramebuffer(fb) : -1;
}

    static int
    reset(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->resetFramebuffer(fb) : -1;
}

    static int
    viewport(imgstream_fb_t *fb, int left, int top, int right, int bottom,
	     void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferViewport(fb, left, top, right, bottom) : -1;
}

    static int
    view(imgstream_fb_t *fb, const imgstream_fb_view_t *view, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferView(fb, view) : -1;
}

    static int
    cursor(imgstream_fb_t *fb, const imgstream_fb_cursor_t *cursor, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferCursor(fb, cursor) : -1;
}

    static int
    scursor(imgstream_fb_t *fb, int mode, int x, int y, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferScreenCursor(fb, mode, x, y) : -1;
}

    static int
    setcursor(imgstream_fb_t *fb, const unsigned char *bits, int xbits, int ybits,
	      int xorig, int yorig, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferCursorShape(fb, bits, xbits, ybits,
	    xorig, yorig) : -1;
}

    static int
    poll(imgstream_fb_t *UNUSED(fb), void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->poll(NULL) : -1;
}

    static long
    pollRate(const imgstream_fb_t *UNUSED(fb), void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->pollRate() : 0;
}
};

imgstream_fb_t *
bobol_window_host_open_display_framebuffer(BObolWindowHost *host,
	const char *spec, size_t width, size_t height)
{
    if (!host)
	return NULL;

    struct imgstream_fb_display_host displayHost;
    memset(&displayHost, 0, sizeof(displayHost));
    displayHost.open = BObolFramebufferBridge::open;
    displayHost.close = BObolFramebufferBridge::close;
    displayHost.flush = BObolFramebufferBridge::flush;
    displayHost.reset = BObolFramebufferBridge::reset;
    displayHost.viewport = BObolFramebufferBridge::viewport;
    displayHost.view = BObolFramebufferBridge::view;
    displayHost.cursor = BObolFramebufferBridge::cursor;
    displayHost.scursor = BObolFramebufferBridge::scursor;
    displayHost.setcursor = BObolFramebufferBridge::setcursor;
    displayHost.poll = BObolFramebufferBridge::poll;
    displayHost.poll_rate = BObolFramebufferBridge::pollRate;
    return imgstream_fb_open_display(spec, width, height, &displayHost, host);
}
